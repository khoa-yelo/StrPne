"""
Build a Genome-Scale Metabolic Model (GEM) from Excel data using CobraPy.

This script reads reaction and metabolite data from an Excel file and constructs
a complete COBRA model with reactions, GPRs, metabolites, and medium definition.
"""

import pandas as pd
import re
from cobra import Model, Reaction, Metabolite
from cobra.io import save_json_model, write_sbml_model


def parse_reaction_equation(equation_str):
    """
    Parse a reaction equation string into reactants and products.

    Handles formats like:
    - "atp[c] + h2o[c] -> adp[c] + pi[c] + h[c]"
    - "g1p[c] <=> g6p[c]"
    - "glc-D[e] <=>  " (exchange reactions)
    - "9.6e-03 d12dg-LLA[c]" (scientific notation coefficients)

    Returns: dict mapping metabolite_id to coefficient (negative for reactants)
    """
    if not equation_str or pd.isna(equation_str):
        return {}

    equation_str = str(equation_str).strip()

    # Determine reversibility and split
    if "<==>" in equation_str:
        left, right = equation_str.split("<==>")
    elif "<=>" in equation_str:
        left, right = equation_str.split("<=>")
    elif "<->" in equation_str:
        left, right = equation_str.split("<->")
    elif "->" in equation_str:
        left, right = equation_str.split("->")
    elif "-->" in equation_str:
        left, right = equation_str.split("-->")
    else:
        # Single metabolite (exchange reaction without arrow)
        return {}

    metabolites = {}

    def parse_side(side_str, sign):
        """Parse one side of the equation."""
        result = {}
        if not side_str or side_str.strip() == "":
            return result

        side_str = side_str.strip()
        if not side_str:
            return result

        # Split by + (careful with scientific notation)
        components = re.split(r"\s*\+\s*", side_str)

        for comp in components:
            comp = comp.strip()
            if not comp:
                continue

            # IMPORTANT: Only treat leading numbers as coefficients if followed by WHITESPACE
            # This prevents "10fthf[c]" from being split as coefficient 10 + metabolite "fthf[c]"
            # Metabolite names can start with numbers (e.g., 10fthf, 5mthf, 2pg, 3pg)

            # Pattern: coefficient (int, float, or scientific) followed by WHITESPACE and metabolite
            # The \s+ (one or more whitespace) is critical to distinguish coefficient from metabolite name
            sci_match = re.match(r"^(\d+\.?\d*(?:[eE][+-]?\d+)?)\s+(.+)$", comp)
            if sci_match:
                coef = float(sci_match.group(1))
                met_id = sci_match.group(2).strip()
            else:
                # No coefficient found (or number is part of metabolite name)
                # Treat entire component as metabolite with coefficient 1
                coef = 1.0
                met_id = comp.strip()

            if met_id:
                result[met_id] = sign * coef

        return result

    # Parse reactants (negative) and products (positive)
    reactants = parse_side(left, -1)
    products = parse_side(right, 1)

    # Combine
    metabolites.update(reactants)
    for met_id, coef in products.items():
        if met_id in metabolites:
            metabolites[met_id] += coef
        else:
            metabolites[met_id] = coef

    return metabolites


def get_compartment_from_id(met_id):
    """Extract compartment from metabolite ID like 'atp[c]' -> 'c'"""
    match = re.search(r"\[(\w+)\]$", met_id)
    if match:
        return match.group(1)
    return "c"  # Default to cytosol


def clean_metabolite_id(met_id):
    """Clean metabolite ID for COBRA compatibility."""
    # Replace problematic characters
    met_id = met_id.replace(" ", "_")
    met_id = met_id.replace("-", "_")
    met_id = met_id.replace("[", "_")
    met_id = met_id.replace("]", "")
    met_id = met_id.replace("(", "_")
    met_id = met_id.replace(")", "")
    return met_id


def build_model_from_excel(excel_file, model_id="D39_model"):
    """
    Build a COBRA model from Excel file containing Reactions and Metabolites sheets.

    Parameters:
    -----------
    excel_file : str
        Path to Excel file with 'Reactions' and 'Metabolites' sheets
    model_id : str
        ID for the model

    Returns:
    --------
    cobra.Model : The constructed model
    """

    print(f"Reading Excel file: {excel_file}")

    # Read Excel sheets
    reactions_df = pd.read_excel(excel_file, sheet_name="Reactions")
    metabolites_df = pd.read_excel(excel_file, sheet_name="Metabolites")

    print(f"Found {len(reactions_df)} reactions and {len(metabolites_df)} metabolites")

    # Create model
    model = Model(model_id)
    model.name = "Streptococcus pneumoniae D39 Genome-Scale Metabolic Model"

    # -------------------------------------------------------------------------
    # Step 1: Create metabolites dictionary from Metabolites sheet
    # -------------------------------------------------------------------------
    print("Creating metabolites...")
    metabolite_info = {}

    for _, row in metabolites_df.iterrows():
        met_abbr = str(row["Abbreviation"]).strip()
        if pd.isna(met_abbr) or not met_abbr:
            continue

        metabolite_info[met_abbr] = {
            "name": str(row["Name"]) if not pd.isna(row["Name"]) else met_abbr,
            "formula": (
                str(row["Formula (charged)"])
                if not pd.isna(row["Formula (charged)"])
                else ""
            ),
            "charge": int(row["Charge"]) if not pd.isna(row["Charge"]) else 0,
            "compartment": get_compartment_from_id(met_abbr),
        }

        # Add annotations
        annotations = {}
        if not pd.isna(row.get("KEGG ID", "")):
            annotations["kegg.compound"] = [str(row["KEGG ID"])]
        if not pd.isna(row.get("ChEBI ID", "")):
            annotations["chebi"] = [str(row["ChEBI ID"])]
        if not pd.isna(row.get("PubChem ID", "")):
            annotations["pubchem.compound"] = [str(row["PubChem ID"])]
        metabolite_info[met_abbr]["annotation"] = annotations

    # -------------------------------------------------------------------------
    # Step 2: Create all metabolites (including ones found in reactions)
    # -------------------------------------------------------------------------
    metabolites_dict = {}  # Store created Metabolite objects

    def get_or_create_metabolite(met_id):
        """Get existing metabolite or create a new one."""
        if met_id in metabolites_dict:
            return metabolites_dict[met_id]

        # Create new metabolite
        cobra_id = clean_metabolite_id(met_id)
        compartment = get_compartment_from_id(met_id)

        if met_id in metabolite_info:
            info = metabolite_info[met_id]
            met = Metabolite(
                id=cobra_id,
                name=info["name"],
                formula=info["formula"] if info["formula"] != "nan" else "",
                charge=info["charge"],
                compartment=compartment,
            )
            if info["annotation"]:
                met.annotation = info["annotation"]
        else:
            # Metabolite not in metabolites sheet, create with minimal info
            met = Metabolite(id=cobra_id, name=met_id, compartment=compartment)

        metabolites_dict[met_id] = met
        return met

    # -------------------------------------------------------------------------
    # Step 3: Create reactions
    # -------------------------------------------------------------------------
    print("Creating reactions...")
    reactions_added = 0
    exchange_reactions = []
    biomass_reaction = None

    for _, row in reactions_df.iterrows():
        rxn_id = str(row["Abbreviation"]).strip()
        if pd.isna(rxn_id) or not rxn_id:
            continue

        # Clean reaction ID
        rxn_id_clean = rxn_id.replace(" ", "_").replace("-", "_")

        # Create reaction
        rxn = Reaction(rxn_id_clean)
        rxn.name = str(row["Name"]) if not pd.isna(row["Name"]) else rxn_id

        # Set bounds
        lower_bound = (
            float(row["Lower bound"]) if not pd.isna(row["Lower bound"]) else 0.0
        )
        upper_bound = (
            float(row["Upper bound"]) if not pd.isna(row["Upper bound"]) else 1000.0
        )
        rxn.lower_bound = lower_bound
        rxn.upper_bound = upper_bound

        # Parse reaction equation and add metabolites
        equation = str(row["Reaction"]) if not pd.isna(row["Reaction"]) else ""
        metabolite_coeffs = parse_reaction_equation(equation)

        if not metabolite_coeffs:
            # Skip reactions without valid equations
            continue

        # Build metabolites dictionary for the reaction
        rxn_mets = {}
        for met_id, coef in metabolite_coeffs.items():
            met_obj = get_or_create_metabolite(met_id)
            rxn_mets[met_obj] = coef

        try:
            rxn.add_metabolites(rxn_mets)
        except Exception as e:
            print(f"Warning: Could not add metabolites to reaction {rxn_id}: {e}")
            continue

        # Set GPR (Gene-Protein-Reaction rule)
        gpr = str(row["GPR"]) if not pd.isna(row["GPR"]) else ""
        if gpr and gpr != "nan":
            try:
                rxn.gene_reaction_rule = gpr
            except Exception as e:
                print(f"Warning: Invalid GPR for {rxn_id}: {gpr}")

        # Set subsystem
        subsystem = str(row["Subsystem"]) if not pd.isna(row["Subsystem"]) else ""
        if subsystem and subsystem != "nan":
            rxn.subsystem = subsystem.strip()

        # Add annotations
        annotations = {}
        if not pd.isna(row.get("EC. Number", "")):
            ec = str(row["EC. Number"])
            if ec != "nan":
                annotations["ec-code"] = [ec]
        if annotations:
            rxn.annotation = annotations

        # Track exchange reactions and biomass
        if rxn_id.startswith("EX") or rxn_id.startswith("Ex"):
            exchange_reactions.append(rxn)

        if "bio" in rxn_id.lower() and (
            "biomass" in rxn.name.lower() or rxn_id.lower() == "bior"
        ):
            biomass_reaction = rxn

        model.add_reactions([rxn])
        reactions_added += 1

    print(f"Added {reactions_added} reactions to model")
    print(f"Found {len(exchange_reactions)} exchange reactions")

    # -------------------------------------------------------------------------
    # Step 4: Set objective function
    # -------------------------------------------------------------------------
    print("Setting objective function...")

    # Try to find biomass reaction
    if biomass_reaction:
        model.objective = biomass_reaction
        print(f"Objective set to: {biomass_reaction.id}")
    else:
        # Search for biomass reaction by name
        for rxn in model.reactions:
            if "biomass" in rxn.id.lower() or "bior" in rxn.id.lower():
                model.objective = rxn
                print(f"Objective set to: {rxn.id}")
                break
        else:
            print(
                "Warning: No biomass reaction found. You may need to set objective manually."
            )

    # -------------------------------------------------------------------------
    # Step 5: Define medium (exchange reactions with negative lower bounds)
    # -------------------------------------------------------------------------
    print("Configuring medium...")

    medium = {}
    for rxn in model.reactions:
        if rxn.id.startswith("EX") or rxn.id.startswith("Ex"):
            if rxn.lower_bound < 0:
                medium[rxn.id] = abs(rxn.lower_bound)

    if medium:
        try:
            model.medium = medium
            print(f"Set {len(medium)} exchange reactions as medium components")
        except Exception as e:
            print(f"Note: Could not set medium directly: {e}")

    return model


def validate_model(model):
    """Validate the model and print summary statistics."""
    print("\n" + "=" * 60)
    print("MODEL VALIDATION")
    print("=" * 60)

    print(f"\nModel ID: {model.id}")
    print(f"Model Name: {model.name}")
    print(f"\nComponents:")
    print(f"  - Reactions: {len(model.reactions)}")
    print(f"  - Metabolites: {len(model.metabolites)}")
    print(f"  - Genes: {len(model.genes)}")

    # Count reactions by type
    exchange_rxns = [
        r for r in model.reactions if r.id.startswith("EX") or r.id.startswith("Ex")
    ]
    transport_rxns = [
        r
        for r in model.reactions
        if "transport" in r.name.lower() or "t" in r.subsystem.lower()
        if r.subsystem
    ]

    print(f"\nReaction types:")
    print(f"  - Exchange reactions: {len(exchange_rxns)}")
    print(
        f"  - Reactions with GPR: {len([r for r in model.reactions if r.gene_reaction_rule])}"
    )

    # Count compartments
    compartments = set(m.compartment for m in model.metabolites)
    print(f"\nCompartments: {compartments}")

    # Medium analysis
    print(f"\nMedium (uptake reactions):")
    medium_count = 0
    for rxn in model.reactions:
        if (rxn.id.startswith("EX") or rxn.id.startswith("Ex")) and rxn.lower_bound < 0:
            medium_count += 1
    print(f"  - {medium_count} exchange reactions allow uptake")

    # Objective
    print(f"\nObjective:")
    print(f"  - {model.objective.expression}")

    return True


def test_model_optimization(model):
    """Test model optimization (FBA)."""
    print("\n" + "=" * 60)
    print("FLUX BALANCE ANALYSIS (FBA)")
    print("=" * 60)

    try:
        solution = model.optimize()

        print(f"\nOptimization status: {solution.status}")
        print(f"Objective value: {solution.objective_value:.6f}")

        if solution.status == "optimal" and solution.objective_value > 0:
            print("\n✓ Model optimizes successfully!")

            # Show some active fluxes
            print("\nTop 10 reactions by absolute flux:")
            fluxes = solution.fluxes.abs().sort_values(ascending=False)
            for rxn_id in fluxes.head(10).index:
                rxn = model.reactions.get_by_id(rxn_id)
                print(f"  {rxn_id}: {solution.fluxes[rxn_id]:.4f} ({rxn.name[:50]}...)")

            # Show exchange fluxes
            print("\nActive exchange reactions:")
            exchange_count = 0
            for rxn in model.reactions:
                if rxn.id.startswith("EX") or rxn.id.startswith("Ex"):
                    flux = solution.fluxes[rxn.id]
                    if abs(flux) > 1e-6:
                        exchange_count += 1
                        if exchange_count <= 15:
                            direction = "uptake" if flux < 0 else "secretion"
                            print(f"  {rxn.id}: {flux:.4f} ({direction})")
            if exchange_count > 15:
                print(f"  ... and {exchange_count - 15} more")
        else:
            print("\n✗ Model did not find optimal solution")
            print(
                "  This may be due to missing exchange reactions or medium constraints"
            )

        return solution

    except Exception as e:
        print(f"\n✗ Optimization failed: {e}")
        return None


def main():
    """Main function to build and test the model."""

    # Build model
    excel_file = "Additional_file_1.xlsx"
    model = build_model_from_excel(excel_file, model_id="iSPD39")

    # Validate
    validate_model(model)

    # Test optimization
    solution = test_model_optimization(model)

    # Save model
    print("\n" + "=" * 60)
    print("SAVING MODEL")
    print("=" * 60)

    # Save as JSON
    json_file = "iSPD39_model.json"
    save_json_model(model, json_file)
    print(f"Model saved to: {json_file}")

    # Save as SBML
    sbml_file = "iSPD39_model.xml"
    write_sbml_model(model, sbml_file)
    print(f"Model saved to: {sbml_file}")

    print("\n" + "=" * 60)
    print("DONE!")
    print("=" * 60)

    return model, solution


if __name__ == "__main__":
    model, solution = main()
