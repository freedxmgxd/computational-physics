from ast import Str
from mp_api.client import MPRester
from itertools import product

API_KEY = "bMpDCd5C2QI8845pHRTxb4JMk496t9Ru"

A_elements = ["Ca", "Sr", "Ba", "K", "Na", "Pb"]
B_elements = ["Ti", "Zr", "Nb", "Ta", "Mn", "Fe"]
X_element = "O"  # Oxigênio como o elemento comum X

Band_gaps = []
Structures = []

band_gap_min = 1.0
band_gap_max = 3.0


def query_perovskites_band_gap():
    with MPRester(API_KEY) as mpr:
        try:
            chemsys_filter = [
                f"{a}-{b}-{X_element}" for a, b in product(A_elements, B_elements)
            ]
            perovskites = mpr.materials.summary.search(
                chemsys=chemsys_filter,
                fields=["material_id", "formula_pretty", "band_gap", "structure"],
            )

            for material in perovskites:
                if material.band_gap < band_gap_min: # or material.band_gap > band_gap_max:
                    continue

                print(f"Material ID: {material.material_id}")
                print(f"Formula: {material.formula_pretty}")
                print(f"Band Gap: {material.band_gap} eV")
                print("-" * 40)

                Band_gaps.append(material.band_gap)
                # print(f"Structure: {material.structure}")
                Structures.append(material.structure)

            print(f"Found {len(Band_gaps)} perovskites with band gap data:")
        except Exception as e:
            print(f"Error during query: {e}")


query_perovskites_band_gap()
