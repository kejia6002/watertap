#  this is a file to practice pressure drop methods in ED
#  import modules

from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    TerminationCondition,
    value,
    Constraint,
    Var,
    Objective,
    Expression,
    assert_optimal_termination,
    log,
)
from pyomo.environ import units as pyunits
from pyomo.util.check_units import (
    assert_units_consistent,
    assert_units_equivalent,
    check_units_equivalent,
)
import pyomo.util.infeasible as infeas
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_unfixed_variables,
    number_total_constraints,
    report_statistics,
    variables_near_bounds_generator,
    total_constraints_set,
)
from idaes.core.util import DiagnosticsToolbox
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
import idaes.core.util.model_diagnostics as m_diag
import idaes.logger as idaeslog

from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models.electrodialysis_1D import (
    ElectricalOperationMode,
    Electrodialysis1D,
    PressureDropMethod,
    FrictionFactorMethod,
    HydraulicDiameterMethod,
)


def main():
    solver = get_solver()
    # create model, flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # create dict to define ions (the prop pack of Adam requires this)
    ion_dict = {
        "solute_list": ["Na_+", "Cl_-"],
        "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3},
        "elec_mobility_data": {("Liq", "Na_+"): 5.19e-8, ("Liq", "Cl_-"): 7.92e-8},
        "charge": {"Na_+": 1, "Cl_-": -1},
    }
    m.obj = Objective(expr=0)

    # attach prop pack to flowsheet
    m.fs.properties = MCASParameterBlock(**ion_dict)
    # build the unit model, pass prop pack to the model
    # m.fs.unit = Electrodialysis1D(default = {"property_package": m.fs.properties,
    # "operation_mode": "Constant_Voltage","finite_elements": 20})
    m.fs.EDstack = Electrodialysis1D(
        property_package=m.fs.properties,
        operation_mode=ElectricalOperationMode.Constant_Voltage,
        has_pressure_change=True,
        finite_elements=10,
        has_nonohmic_potential_membrane=False,
        has_Nernst_diffusion_layer=False,
        # pressure_drop_method=PressureDropMethod.Darcy_Weisbach,
        pressure_drop_method=PressureDropMethod.experimental,
        friction_factor_method=FrictionFactorMethod.Gurreri,
        hydraulic_diameter_method=HydraulicDiameterMethod.fixed,
    )
    assert_units_consistent(m)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 0.1, index=("Liq", "H20")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 5e2, index=("Liq", "Na_+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 5e2, index=("Liq", "Cl_-")
    )
    iscale.set_scaling_factor(m.fs.EDstack.cell_width, 5)
    iscale.set_scaling_factor(m.fs.EDstack.cell_length, 5)
    iscale.set_scaling_factor(m.fs.EDstack.cell_pair_num, 5)
    iscale.set_scaling_factor(m.fs.EDstack.voltage_applied, 5)

    iscale.calculate_scaling_factors(m.fs.EDstack)

    # specify the model by fixing parameters
    # fix state variables in diluate and concentrate
    # 3 flow_mol_phase_comp for each chamber, 6

    m.fs.EDstack.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(110)
    m.fs.EDstack.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0.342)
    m.fs.EDstack.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(0.342)
    m.fs.EDstack.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(110)
    m.fs.EDstack.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0.342)
    m.fs.EDstack.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(0.342)

    # 2 temperature and pressure for each chamber ,4
    m.fs.EDstack.inlet_diluate.pressure[0].fix(5003485)
    m.fs.EDstack.inlet_diluate.temperature.fix(298.15)
    m.fs.EDstack.inlet_concentrate.pressure[0].fix(5003485)
    m.fs.EDstack.inlet_concentrate.temperature.fix(298.15)

    #  operation condition (voltage or current)
    m.fs.EDstack.voltage_applied.fix(88.5)
    m.fs.EDstack.cell_width.fix(0.197)
    m.fs.EDstack.cell_length.fix(1.69)
    m.fs.EDstack.cell_pair_num.fix(250)

    m.fs.EDstack.current_utilization.fix(1)
    m.fs.EDstack.channel_height.fix(2.7e-4)

    m.fs.EDstack.electrodes_resistance.fix(0)
    m.fs.EDstack.spacer_porosity.fix(0.83)

    m.fs.EDstack.water_trans_number_membrane["cem"].fix(5.8)
    m.fs.EDstack.water_trans_number_membrane["aem"].fix(4.3)
    m.fs.EDstack.water_permeability_membrane["cem"].fix(2.16e-14)
    m.fs.EDstack.water_permeability_membrane["aem"].fix(1.75e-14)
    m.fs.EDstack.membrane_areal_resistance["cem"].fix(1.89e-4)
    m.fs.EDstack.membrane_areal_resistance["aem"].fix(1.77e-4)
    m.fs.EDstack.membrane_thickness["aem"].fix(1.3e-4)
    m.fs.EDstack.membrane_thickness["cem"].fix(1.3e-4)
    m.fs.EDstack.solute_diffusivity_membrane["cem", "Na_+"].fix(1.8e-10)
    m.fs.EDstack.solute_diffusivity_membrane["aem", "Na_+"].fix(1.25e-10)
    m.fs.EDstack.solute_diffusivity_membrane["cem", "Cl_-"].fix(1.8e-10)
    m.fs.EDstack.solute_diffusivity_membrane["aem", "Cl_-"].fix(1.25e-10)
    m.fs.EDstack.ion_trans_number_membrane["cem", "Na_+"].fix(1)
    m.fs.EDstack.ion_trans_number_membrane["aem", "Na_+"].fix(0)
    m.fs.EDstack.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
    m.fs.EDstack.ion_trans_number_membrane["aem", "Cl_-"].fix(1)

    # fix for frictional factor calculation
    m.fs.EDstack.diffus_mass.fix(1.6e-9) if hasattr(m.fs.EDstack, "diffus_mass") else 0
    #
    # fix for the hydraulic diameter calculation in spacer_specficed method
    # m.fs.EDstack.spacer_specific_area.fix(10700) if hasattr(
    #     m.fs.EDstack, "spacer_specific_area"
    # ) else 0

    m.fs.EDstack.hydraulic_diameter.fix(4.47e-4) if hasattr(
        m.fs.EDstack, "hydraulic_diameter"
    ) else 0
    # m.fs.EDstack.hydraulic_diameter.fix(3.47e-4) if hasattr(
    #     m.fs.EDstack, "hydraulic_diameter") else 0
    # m.fs.EDstack.hydraulic_diameter.fix(2*7.1e-4)

    print("nnnnnn") if hasattr(m.fs.EDstack, "hydraulic_diameter") else print("yyyy")

    # hydraulic_diameter_conventional = (
    #         2*m.fs.EDstack.channel_height() *
    #         m.fs.EDstack.cell_width() *
    #         m.fs.EDstack.spacer_porosity() *
    #         (m.fs.EDstack.channel_height() +
    #          m.fs.EDstack.cell_width()) ** -1
    # )
    # spacer_specific_area = 1e4
    #
    # hydraulic_diameter_specified = (
    #         4 * m.fs.EDstack.spacer_porosity() *
    #         (2 * m.fs.EDstack.channel_height() ** -1 +
    #         (1 - m.fs.EDstack.spacer_porosity()) *
    #          spacer_specific_area)**-1)
    # friction_factor_Gurreri = (4*50.6*m.fs.EDstack.spacer_porosity() ** -7.06*m.fs.EDstack.N_Re() ** -1)
    # friction_factor_Kuroda = (4 * 9.6 * m.fs.EDstack.spacer_porosity() ** -1 * m.fs.EDstack.N_Re() ** -0.5)
    # print("friction_factor_Gurreri", friction_factor_Gurreri)
    # print("friction_factor_Kuroda", friction_factor_Kuroda)
    # print("hydraulic diameter using specified methods:",
    #       hydraulic_diameter_specified
    #       )
    #
    # print("hydraulic diameter using conventional methods:",
    #       hydraulic_diameter_conventional
    #       )

    # m.fs.EDstack.friction_factor.fix(20)
    # m.fs.EDstack.pressure_drop.fix(1.6e6)

    print("DOF IS ")
    print(degrees_of_freedom(m))
    solver = get_solver()
    m.fs.EDstack.initialize(optarg=solver.options)
    m.fs.EDstack.report()

    results = solver.solve(m)
    print(results.solver.termination_condition)
    return m


if __name__ == "__main__":
    m = main()
