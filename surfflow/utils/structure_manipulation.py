from typing import List, Union

from emmet.core.mpid import MPID
from mp_api.client import MPRester
from pymatgen.core import Structure
from pymatgen.core.surface import Slab
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as SgA

from surfflow.utils.db_tools import VaspDB
from htflow_utils.misc_tools import transfer_average_magmoms


def get_bulk_from_mp(mpid: str, requested_fields: List[str] = None) -> Structure:
    """
    Queries the conventional standard bulk structure for a material given by its
    MaterialsProject ID from the MaterialsProject API.

    :param mpid: Unique MaterialsProject ID of the material.
    :type mpid: str

    :param requested_fields: List of fields to request from the MaterialsProject API.
    :type requested_fields: List[str], optional

    :return: pymatgen.core.structure.Structure
    :rtype: Structure
    """
    with MPRester() as mpr:
        allowed_fields = mpr.summary.available_fields
        if requested_fields is None:
            requested_fields = []
        else:
            requested_fields = [
                field for field in requested_fields if field in allowed_fields
            ]
        requested_fields.extend(["structure", "formula_pretty"])

        doc = mpr.summary.search(material_ids=[MPID(mpid)], fields=requested_fields)[0]
    return doc


def get_conv_bulk_from_db(
    mpid: str, coll: str, db_file: str = "auto", high_level: Union[str, bool] = True
) -> Structure:
    """
    Queries the conventional standard bulk structure for a material given by its
    MaterialsProject ID from the database. If not found, queries it from the
    MaterialsProject API.

    :param mpid: Unique MaterialsProject ID of the material.
    :type mpid: str

    :param coll: Collection in which to search for the conventional standard bulk
        structure in the database.
    :type coll: str

    :param db_file: Full path of the db.json file. Defaults to 'auto', in which case the path
        is checked from the environment variable FW_CONFIG_FILE.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
        Defaults to True, in which case the value in the db.json file is used.
    :type high_level: bool, str, optional

    :raises ValueError: Either the fully optimized conventional standard bulk structure,
        or the one directly loaded from MP should already be in the database.

    :return: pymatgen.core.structure.Structure
    :rtype: Structure
    """
    nav = VaspDB(db_file, high_level)
    bulk_dict = nav.find_data(coll, {"mpid": mpid})

    if bulk_dict:
        try:
            bulk = Structure.from_dict(bulk_dict["structure_fromMP_conv"])
        except KeyError:
            try:
                bulk = Structure.from_dict(bulk_dict["structure_fromMP"])
            except KeyError:
                raise ValueError(
                    "Either the fully optimized primitive bulk structure, "
                    "or the one directly loaded from MP should already be "
                    "in the database."
                )
            else:
                bulk = SgA(bulk).get_conventional_standard_structure(
                    keep_site_properties=True
                )
    else:
        add_bulk_to_db(mpid, coll, db_file, high_level)
        bulk = get_conv_bulk_from_db(mpid, coll, db_file, high_level)

    return bulk


def add_bulk_to_db(
    mpid: str,
    coll: str,
    db_file: str = "auto",
    high_level: Union[str, bool] = True,
    custom_data: dict = None,
) -> None:
    """
    Adds the conventional standard bulk structure for a material given by its
    MaterialsProject ID to the database.

    :param custom_data:
    :param mpid: Unique MaterialsProject ID of the material.
    :type mpid: str

    :param coll: Collection in which to search for the conventional standard bulk
        structure in the database.
    :type coll: str

    :param db_file: Full path of the db.json file. Defaults to 'auto', in which case the path
        is checked from the environment variable FW_CONFIG_FILE.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
        Defaults to True, in which case the value in the db.json file is used.
    :type high_level: bool, str, optional

    :return: None
    """
    nav = VaspDB(db_file=db_file, high_level=high_level)

    bulk_data = nav.find_data(coll, {"mpid": mpid})
    if bulk_data:
        print(f"Bulk structure for {mpid} already in database.")
        return

    bulk_doc = get_bulk_from_mp(
        mpid=mpid, requested_fields=["band_gap", "is_magnetic", "ordering"]
    )
    bulk_struct = bulk_doc.structure
    prim_bulk_struct = bulk_struct.get_primitive_structure(use_site_props=True)
    conv_bulk_struct = SgA(prim_bulk_struct).get_conventional_standard_structure(keep_site_properties=True)

    conv_bulk_struct = transfer_average_magmoms(bulk_struct, conv_bulk_struct)

    data = {
        "mpid": mpid,
        "formula": bulk_doc.formula_pretty,
        "structure_fromMP": bulk_struct.as_dict(),
        "structure_fromMP_conv": conv_bulk_struct.as_dict(),
        "primitive_structure": prim_bulk_struct.as_dict(),
        "band_gap": bulk_doc.band_gap,
        "is_magnetic": bulk_doc.is_magnetic,
        "ordering": bulk_doc.ordering,
    }

    if custom_data:
        data.update(custom_data)

    if comp_params := data.get("comp_parameters"):
        comp_params.update({"is_metal": bulk_doc.band_gap < 0.3})

    nav.update_data(
        collection=coll,
        fltr={"mpid": mpid},
        new_values={"$set": data},
        upsert=True,
    )


def find_slabs(
    fltr: dict,
    collection: str,
    db_file: str = "auto",
    high_level: Union[str, bool] = True,
):
    """
    Return a slab and a string showing if it is fully relaxed or just pre-relaxed.

    This function is useful to fully relax slabs after the surface energy is
    calculated. If the slab is asymmetric, top and bottom parts are relaxed
    separately. They can be combined to form a pre-relaxed slab.
    For symmetrical slabs the surface energy workflow already fully relaxes
    them, thus it can just be grabbed from there.

    :param fltr: Filter to find the slab in the database.
    :type fltr: dict

    :param collection: Collection in which to search for the slab in the database.
    :type collection: str

    :param db_file: Full path of the db.json file. Defaults to 'auto', in which case the path
        is checked from the environment variable FW_CONFIG_FILE.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
        Defaults to True, in which case the value in the db.json file is used.
    :type high_level: bool, str, optional

    :raises ValueError: If the slab is not found in the database.

    :return: Slab structure and a string showing if it is fully relaxed or just pre-relaxed.
    :rtype: str, pymatgen.core.structure.Structure
    """
    # find the calculation dictionary for the correct material, miller index,
    # and termination.
    nav = VaspDB(db_file=db_file, high_level=high_level)
    data_dict = nav.find_data(collection=collection, fltr=fltr)
    calcs_dict = data_dict["calcs"]

    # Check if the fully relaxed slab is already there and return it if true.
    try:
        fully_relaxed_slab = Structure.from_dict(
            calcs_dict["slab_relax"]["output"]["structure"]
        )
        return "fully_relaxed", fully_relaxed_slab
    except KeyError:
        pass

    # Try to find partially relaxed slabs, and combine them using the
    # get_prerelaxed_slab function.
    try:
        relaxed_top = Structure.from_dict(
            calcs_dict["slab_top_fixed_relax"]["output"]["structure"]
        )
        relaxed_bot = Structure.from_dict(
            calcs_dict["slab_bot_fixed_relax"]["output"]["structure"]
        )
        unrelaxed_slab = Slab.from_dict(calcs_dict["slab_static"]["structure"])
        pre_relaxed_slab = get_prerelaxed_slab(relaxed_top, relaxed_bot, unrelaxed_slab)
        return "not_fully_relaxed", pre_relaxed_slab
    except KeyError:
        pass

    # If neither can be found, return None.
    return "nothing_found", None


def get_prerelaxed_slab(
    relaxed_top: Structure, relaxed_bot: Structure, unrelaxed_slab: Slab
) -> Slab:
    """
    Combine partially relaxed slabs to form one pre-relaxed slab as a good starting guess.

    This function is useful to fully relax slabs after the surface energy is
    calculated. If the slab is asymmetric, top and bottom parts are relaxed
    separately. They can be combined to form a pre-relaxed slab. For symmetrical
    slabs the surface energy workflow already fully relaxes them, thus it can
    just be grabbed from there.

    :param relaxed_top: Slab whose top sites are relaxed.
    :type relaxed_top: pymatgen Structure, or Slab

    :param relaxed_bot: Slab whose bottom sites are relaxed.
    :type relaxed_bot: pymatgen Structure, or Slab

    :param unrelaxed_slab: unrelaxed slab used to get site properties and
        potentially sites that have not been relaxed in either slab.
    :type unrelaxed_slab: pymatgen Slab

    :return: Slab with relaxed sites on top and bottom.
    :rtype: pymatgen Slab
    """
    # copy the unrelaxed slab and remove the copy's sites
    new_slab = unrelaxed_slab.copy()
    new_slab.sites.clear()

    # assign relaxed sites and assign properties from original one.
    # sites are assumed to be relaxed if their selective dynamics flags are
    # all true
    for i, s in enumerate(unrelaxed_slab.sites):
        top_site = relaxed_top.sites[i]
        bot_site = relaxed_bot.sites[i]
        if all(top_site.properties["selective_dynamics"]):
            new_slab.append(
                species=top_site.species,
                coords=top_site.frac_coords,
                properties=s.properties,
            )
        elif all(bot_site.properties["selective_dynamics"]):
            new_slab.append(
                species=bot_site.species,
                coords=bot_site.frac_coords,
                properties=s.properties,
            )
        # if we have no relaxed sites in either slab, take the unrelaxed ones.
        else:
            new_slab.append(
                species=s.species, coords=s.frac_coords, properties=s.properties
            )

    return new_slab
