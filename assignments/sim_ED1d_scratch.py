#  this is a file to sim ED_1d model with only an ED stack
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
from watertap.costing.watertap_costing_package import WaterTAPCosting


def main():

    #  create a concrete model from pyomo
    m = ConcreteModel()
    # add flow sheet to this model from idaes
    m.fs = FlowsheetBlock(dynamic=False)
    ion_dict = {
        "solute_list": ["Na_+", "Cl_-"],
        "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3},
        "elec_mobility_data": {("Liq", "Na_+"): 5.19e-8, ("Liq", "Cl_-"): 7.92e-8},
        "charge": {"Na_+": 1, "Cl_-": -1},
    }

    # add properties to this model, in this case it is a multi component aqueous
    # solution property based on given ion_dict
    m.fs.properties = MCASParameterBlock(**ion_dict)
    m.fs.costing = WaterTAPCosting()
    # m.fs.prod = Product(property_package=m.fs.properties)
    # add unit model from watertap
    m.fs.EDstack = Electrodialysis1D(
        property_package=m.fs.properties,
        operation_mode=ElectricalOperationMode.Constant_Voltage,
        has_pressure_change=True,
        finite_elements=10,
        has_nonohmic_potential_membrane=False,
        has_Nernst_diffusion_layer=False,
        pressure_drop_method=PressureDropMethod.Darcy_Weisbach,
        friction_factor_method=FrictionFactorMethod.Gurreri,
        hydraulic_diameter_method=HydraulicDiameterMethod.spacer_specific_area_known,
    )
    m.fs.EDstack.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        # costing_method_argumengts={"cost_electricity_flow": True, "has_rectifier": True},
    )
    # this model does not require connectivity

    # assert the units to see if they are consistent
    assert_units_consistent(m)
    # add variables and constraints: no extra variables, just use the default ones from ED 1d model
    # set the scaling factors for state variables
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 0.1, index=("Liq", "H20")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 5e2, index=("Liq", "Na_+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 5e2, index=("Liq", "Cl_-")
    )
    # call the calculate_scaling_factors function
    # iscale.set_scaling_factor(m.fs.EDstack.cell_width, 5)
    # iscale.set_scaling_factor(m.fs.EDstack.cell_length, 0.1)
    # iscale.set_scaling_factor(m.fs.EDstack.cell_pair_num, 5)
    # iscale.set_scaling_factor(m.fs.EDstack.water_permeability_membrane["cem"], 1e20)
    # iscale.set_scaling_factor(m.fs.EDstack.water_permeability_membrane["aem"], 1e20)
    # for key in m.fs.EDstack.current_density_x.keys():
    # iscale.set_scaling_factor(m.fs.EDstack.current_density_x[key], 1e2)
    # iscale.set_scaling_factor(m.fs.EDstack.voltage_applied_x[key], 1e-6)

    iscale.set_scaling_factor(m.fs.EDstack.cell_width, 1)
    iscale.set_scaling_factor(m.fs.EDstack.cell_length, 1)
    iscale.set_scaling_factor(m.fs.EDstack.cell_pair_num, 0.01)
    iscale.set_scaling_factor(m.fs.EDstack.voltage_applied, 0.2)

    iscale.calculate_scaling_factors(m.fs.EDstack)

    # specify the model by fixing parameters
    # fix state variables in diluate and concentrate
    # 3 flow_mol_phase_comp for each chamber, 6
    m.fs.EDstack.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(28.74)
    m.fs.EDstack.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(4.449e-2)
    m.fs.EDstack.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(4.449e-2)
    m.fs.EDstack.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(28.74)
    m.fs.EDstack.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(4.449e-2)
    m.fs.EDstack.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(4.449e-2)

    # m.fs.EDstack.outlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(28.74)
    m.fs.EDstack.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(2.35e-2)
    # m.fs.EDstack.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(4.449e-3)
    # m.fs.EDstack.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(28.74)
    # m.fs.EDstack.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(5.449e-2)
    # m.fs.EDstack.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(4.449e-3)

    # Corresponding to C_feed = 5g/L, flow rate = 5.2e-4 m3/s
    # the flow rate were lower
    # m.fs.EDstack.diluate.properties[0, 0].flow_vol_phase["Liq"]
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(
        m.fs.EDstack.diluate.properties[0, 0].flow_vol_phase["Liq"]
    )
    m.fs.costing.add_LCOW(m.fs.EDstack.diluate.properties[0, 0].flow_vol_phase["Liq"])
    m.fs.costing.add_specific_energy_consumption(
        m.fs.EDstack.diluate.properties[0, 0].flow_vol_phase["Liq"]
    )
    # 2 temperature and pressure for each chamber ,4
    m.fs.EDstack.inlet_diluate.pressure[0].fix(5003485)
    m.fs.EDstack.inlet_diluate.temperature.fix(298.15)
    m.fs.EDstack.inlet_concentrate.pressure[0].fix(5003485)
    m.fs.EDstack.inlet_concentrate.temperature.fix(298.15)

    #  operation condition (voltage or current)
    m.fs.EDstack.voltage_applied.fix(15)
    # m.fs.EDstack.cell_width.fix(0.197)
    m.fs.EDstack.cell_length.fix(0.69)
    m.fs.EDstack.cell_pair_num.fix(50)

    m.fs.EDstack.current_utilization.fix(1)
    m.fs.EDstack.channel_height.fix(2.7e-4)

    m.fs.EDstack.electrodes_resistance.fix(0)
    m.fs.EDstack.spacer_porosity.fix(0.83)

    # 16
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

    # pressure parameters(4) : pressure drop, diffus_mass, spacer_specific_area, friction_factor
    if m.fs.EDstack.config.pressure_drop_method == PressureDropMethod.experimental:
        m.fs.EDstack.pressure_drop.fix(40000)

    m.fs.EDstack.diffus_mass.fix(1.6e-9) if hasattr(m.fs.EDstack, "diffus_mass") else 0

    m.fs.EDstack.spacer_specific_area.fix(10700) if hasattr(
        m.fs.EDstack, "spacer_specific_area"
    ) else 0
    m.fs.EDstack.friction_factor.fix(
        20
    ) if m.fs.EDstack.config.friction_factor_method == FrictionFactorMethod.fixed else 0

    print("----------------------------------------------")
    print("report model statistics after specifying", report_statistics(m.fs))
    # m.fs.EDstack.pprint()
    print("degree of freedoms", degrees_of_freedom(m))
    # initialize

    print("BADLY SCALED VARS & CONSTRAINTS")
    print("------------------------------------------------------")

    badly_scaled_var_values = {
        var.name: val
        for (var, val) in iscale.badly_scaled_var_generator(
            m, large=100, small=0.01, zero=1e-10
        )
    }
    for j, k in badly_scaled_var_values.items():
        print(j, ":", k)
    print("end of badly scaled var and constraints")
    print("------------------------------------------------------")

    # dt = DiagnosticsToolbox(m)
    # dt.report_structural_issues()
    # # dt.report_numerical_issues()

    solver = get_solver()
    m.fs.EDstack.initialize(optarg=solver.options)
    m.fs.EDstack.report()
    # solve the model
    results = solver.solve(m)
    print(results.solver.termination_condition)
    # optimize the model
    """
    7.1) Unfix some degrees of freedom to provide the problem with decision variables, variable_name.unfix().
    7.2) Add bounds to variables and inequality constraints to constrain solution space, 
    variable_name.setlb(value) and var_name.setub(value)
    7.3) Call a solver and check the termination conditions, see step 6 Solving the Model.
    """
    # m.obj=Objective(expr=m.fs.EDstack.costing.LCOW)
    # m.fs.objective = Objective(expr=m.fs.costing.LCOW)
    # m.fs.EDstack.cell_width.unfix()
    # m.fs.EDstack.cell_length.unfix()
    # m.fs.EDstack.cell_pair_num.unfix()
    # m.fs.EDstack.voltage_applied[0].unfix()git
    # m.fs.EDstack.voltage_applied[0].setlb(0.1)
    # m.fs.EDstack.voltage_applied[0].setub(85)
    # m.fs.EDstack.cell_pair_num.unfix()
    # m.fs.EDstack.cell_pair_num.setlb(1)

    # results = solver.solve(m, tee=False)

    # assert results.solver.termination_condition == TerminationCondition.optimal
    # print(results.solver.termination_condition)

    # m.fs.EDstack.display()

    m.fs.EDstack.cell_width.pprint()
    m.fs.EDstack.cell_length.pprint()
    m.fs.EDstack.cell_pair_num.pprint()
    m.fs.EDstack.voltage_applied.pprint()
    # print("calculator:",(4.49e-2-2.62e-3)/4.49e-2)
    return m


if __name__ == "__main__":
    m = main()
