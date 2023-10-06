from math import pi
from xmlrpc.client import FastParser
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
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_unfixed_variables,
    number_total_constraints,
    report_statistics,
    variables_near_bounds_generator,
    total_constraints_set,
)
import idaes.core.util.model_statistics as stats
from idaes.core.util.constants import Constants
import idaes.core.util.scaling as iscale
import idaes.core.util.model_diagnostics as m_diag
import idaes.logger as idaeslog
import openpyxl
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock

# from ion_DSPMDE_prop_pack_inwork import DSPMDEParameterBlock # testing with the updated property package

from watertap.unit_models.electrodialysis_1D import (
    ElectricalOperationMode,
    Electrodialysis1D,
    PressureDropMethod,
    FrictionFactorMethod,
    HydraulicDiameterMethod,
)

from idaes.core.solvers import get_solver
import pandas as pd


def main():
    solver = get_solver()
    # create model, flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # create dict to define ions (the prop pack of Adam requires this)

    """
    ion_dict = {
        "solute_list": ["Na_+", "Cl_-", "N"],
        "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3, "N": 61.8e-3},
        "elec_mobility_data": {("Liq", "Na_+"): 5.19e-8, ("Liq", "Cl_-"): 7.92e-8},
        "charge": {"Na_+": 1, "Cl_-": -1},
    }
    """
    ion_dict = {
        "solute_list": ["Na_+", "Cl_-"],
        "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3},
        "elec_mobility_data": {("Liq", "Na_+"): 5.19e-8, ("Liq", "Cl_-"): 7.92e-8},
        "charge": {"Na_+": 1, "Cl_-": -1},
        # "visc_d_phase": {"Liq": 890e-6},
    }

    # attach prop pack to flowsheet
    m.fs.properties = MCASParameterBlock(**ion_dict)
    # build the unit model, pass prop pack to the model
    # m.fs.unit = Electrodialysis1D(default = {"property_package": m.fs.properties, "operation_mode": "Constant_Voltage",  "finite_elements": 20})
    m.fs.unit = Electrodialysis1D(
        property_package=m.fs.properties,
        operation_mode=ElectricalOperationMode.Constant_Voltage,
        has_pressure_change=True,
        finite_elements=10,
        has_nonohmic_potential_membrane=False,
        has_Nernst_diffusion_layer=False,
        pressure_drop_method=PressureDropMethod.Darcy_Weisbach,
        friction_factor_method=FrictionFactorMethod.Gurreri,
        hydraulic_diameter_method=HydraulicDiameterMethod.spacer_specific_area_known,
        # limiting_current_density_data= 800,
    )
    m.obj = Objective(expr=0)
    assert_units_consistent(m)

    print("----------------------------------------------")
    # print("report model statistics before specifying", report_statistics(m.fs))
    # m.fs.unit.config.add("hydraulic_diameter_method", HydraulicDiameterMethod.conventional)
    m.fs.unit.inlet_diluate.pressure[0].fix(5003485)
    m.fs.unit.inlet_diluate.temperature.fix(298.15)
    m.fs.unit.inlet_concentrate.pressure[0].fix(5003485)
    m.fs.unit.inlet_concentrate.temperature.fix(298.15)

    # for ind, c in m.fs.unit.diluate.properties.items():
    #     if ind == [0,0]:
    #         c.calculate_state({
    #                     ("flow_vol_phase","Liq"):3.25e-4,
    #                     ("conc_mass_phase_comp",("Liq","Na_+")):3.932,
    #                     ("conc_mass_phase_comp",("Liq","Cl_-")):6.068},
    #                     hold_state=True)
    # for ind, c in m.fs.unit.concentrate.properties.items():
    #     if ind == [0,0]:
    #         c.calculate_state({
    #                     ("flow_vol_phase","Liq"):3.25e-4,
    #                     ("conc_mass_phase_comp",("Liq","Na_+")):3.932,
    #                     ("conc_mass_phase_comp",("Liq","Cl_-")):6.068},
    #                     hold_state=True)

    m.fs.unit.concentrate.properties[0, 0].calculate_state(
        {
            ("flow_vol_phase", "Liq"): 3.25e-4,
            ("conc_mmoass_phase_comp", ("Liq", "Na_+")): 3.932,
            ("conc_mass_phase_comp", ("Liq", "Cl_-")): 6.068,
        },
        hold_state=True,
    )

    # m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(5.777e4)
    # m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(4.449e-8)
    # m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(4.449e-8)
    # m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(5.777e4)
    # m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(4.449e-8)
    # m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(4.449e-8)
    #
    # m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(110)
    # m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0.342)
    # m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(0.342)
    # m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(110)
    # m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0.342)
    # m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(0.342)
    # m.fs.unit.inlet_diluate.flow_vol_phase[0, "Liq", "H2O"].fix(5.2e-4)
    # m.fs.unit.inlet_diluate.conc_mol_phase_comp[0, "Liq", "Na_+"].fix(34.188)
    # m.fs.unit.inlet_diluate.conc_mol_phase_comp[0, "Liq", "Cl_-"].fix(34.188)
    # m.fs.unit.inlet_concentrate.flow_vol_phase[0, "Liq", "H2O"].fix(5.2e-4)
    # m.fs.unit.inlet_concentrate.conc_mol_phase_comp[0, "Liq", "Na_+"].fix(34.188)
    # m.fs.unit.inlet_concentrate.conc_mol_phase_comp[0, "Liq", "Cl_-"].fix(34.188)

    # flow_solute, in molar / s
    # 4.449007529089664e-08
    # flow_water, in molar / s
    # 5777.777777777777

    #
    init_arg = {
        ("flow_vol_phase", ("Liq")): 5.2e-4,
        ("conc_mol_phase_comp", ("Liq", "Na_+")): 34.188,
        ("conc_mol_phase_comp", ("Liq", "Cl_-")): 34.188,
    }  # Corresponding to C_feed = 2g/L
    # m.fs.unit.concentrate.calculate_state(
    #     init_arg,
    #     hold_state=True
    # )
    # m.fs.unit.concentrate.properties[0, 0].calculate_state({
    #     ("flow_vol_phase", "Liq"): 3.25e-4,
    #     ("conc_mass_phase_comp", ("Liq", "Na_+")): 3.932,
    #     ("conc_mass_phase_comp", ("Liq", "Cl_-")): 6.068},
    #     hold_state=True)
    #
    # m.fs.unit.concentrate.properties.calculate_state(
    #     init_arg,
    #     hold_state=True
    # )
    # m.fs.unit.inlet_diluate.properties.calculate_state(
    #     init_arg,
    #     hold_state=True
    # )
    m.fs.unit.water_trans_number_membrane["cem"].fix(5.8)
    m.fs.unit.water_trans_number_membrane["aem"].fix(4.3)
    m.fs.unit.water_permeability_membrane["cem"].fix(2.16e-14)
    m.fs.unit.water_permeability_membrane["aem"].fix(1.75e-14)
    m.fs.unit.voltage_applied.fix(88.5)
    # m.fs.unit.current_applied.fix(8)
    m.fs.unit.electrodes_resistance.fix(0)
    m.fs.unit.cell_pair_num.fix(250)
    m.fs.unit.current_utilization.fix(1)
    m.fs.unit.channel_height.fix(2.7e-4)
    m.fs.unit.membrane_areal_resistance["cem"].fix(1.89e-4)
    m.fs.unit.membrane_areal_resistance["aem"].fix(1.77e-4)
    m.fs.unit.cell_width.fix(0.197)
    m.fs.unit.cell_length.fix(1.68)
    m.fs.unit.membrane_thickness["aem"].fix(1.3e-4)
    m.fs.unit.membrane_thickness["cem"].fix(1.3e-4)
    m.fs.unit.solute_diffusivity_membrane["cem", "Na_+"].fix(1.8e-10)
    m.fs.unit.solute_diffusivity_membrane["aem", "Na_+"].fix(1.25e-10)
    m.fs.unit.solute_diffusivity_membrane["cem", "Cl_-"].fix(1.8e-10)
    m.fs.unit.solute_diffusivity_membrane["aem", "Cl_-"].fix(1.25e-10)
    # m.fs.unit.solute_diffusivity_membrane["cem", "N"].fix(1.8e-10)
    # m.fs.unit.solute_diffusivity_membrane["aem", "N"].fix(1.25e-10)
    m.fs.unit.ion_trans_number_membrane["cem", "Na_+"].fix(1)
    m.fs.unit.ion_trans_number_membrane["aem", "Na_+"].fix(0)
    m.fs.unit.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
    m.fs.unit.ion_trans_number_membrane["aem", "Cl_-"].fix(1)

    if m.fs.unit.config.pressure_drop_method == PressureDropMethod.experimental:
        m.fs.unit.pressure_drop.fix(40000)
    m.fs.unit.spacer_porosity.fix(0.83)

    m.fs.unit.diffus_mass.fix(1.6e-9) if hasattr(m.fs.unit, "diffus_mass") else 0
    # m.fs.unit.hydraulic_diameter.fix(1.5e-3) if hasattr(m.fs.unit, "diffus_mass") else 0
    # m.fs.unit.visc_d.pprint()
    m.fs.unit.spacer_specific_area.fix(10700) if hasattr(
        m.fs.unit, "spacer_specific_area"
    ) else 0
    m.fs.unit.friction_factor.fix(
        20
    ) if m.fs.unit.config.friction_factor_method == FrictionFactorMethod.fixed else 0
    # m.fs.unit.hydraulic_diameter.fix(2*7.1e-4)

    # m.fs.unit.config.property_package.visc_d_phase["Liq"].fix(890e-6)

    print("----------------------------------------------")
    print("report model statistics after specifying", report_statistics(m.fs))

    # m.fs.unit.pprint()
    print("degrees of freedom", degrees_of_freedom(m.fs))
    assert degrees_of_freedom(m.fs) == 0

    # set scaling factors for state vars and call the 'calculate_scaling_factors' function
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 0.1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "Na_+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "Cl_-")
    )
    # m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1e5, index=("Liq", "N"))
    iscale.set_scaling_factor(m.fs.unit.cell_width, 5)
    iscale.set_scaling_factor(m.fs.unit.cell_length, 1)
    iscale.set_scaling_factor(m.fs.unit.cell_pair_num, 0.01)
    iscale.set_scaling_factor(m.fs.unit.voltage_applied, 0.1)
    # m.fs.unit.diluate.display()
    iscale.calculate_scaling_factors(m.fs.unit)
    m.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)
    # m.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

    print("BADLY SCALED VARS & CONSTRAINS")
    badly_scaled_var_values = {
        var.name: val
        for (var, val) in iscale.badly_scaled_var_generator(
            m, large=100, small=0.01, zero=1e-10
        )  # ,
    }
    for j, k in badly_scaled_var_values.items():
        print(j, ":", k)

    results = solver.solve(m, tee=False, symbolic_solver_labels=True)
    print(
        "Solver status: {}, Termination condition: {}".format(
            results.solver.status, results.solver.termination_condition
        )
    )
    print("NEAR-BOUND VAR")
    for i in variables_near_bounds_generator(m):
        print(i.name, i.value)
    dh = m_diag.DegeneracyHunter(m)
    dh.check_residuals(tol=1e-5)

    m.fs.unit.report()
    m.fs.unit.current_density_x.pprint()
    # m.fs.unit.diluate.power_electrical_x.pprint()
    # m.fs.unit.specific_power_electrical.pprint()
    # m.fs.unit.total_areal_resistance_x.pprint()
    if hasattr(m.fs.unit, "pressure_drop"):
        m.fs.unit.pressure_drop.pprint()
    else:
        pass
    m.fs.unit.velocity_diluate.pprint()
    # m.fs.unit.velocity_concentrate.pprint()
    m.fs.unit.diluate.deltaP.pprint() if hasattr(m.fs.unit.diluate, "deltaP") else 0
    m.fs.unit.pressure_drop_total.pprint() if hasattr(
        m.fs.unit, "pressure_drop_total"
    ) else 0

    m.fs.unit.N_Re.pprint() if hasattr(m.fs.unit, "N_Re") else 0
    m.fs.unit.diluate.properties[0, 0].pressure.pprint()
    m.fs.unit.diluate.properties[0, 1].pressure.pprint()
    m.fs.unit.concentrate.properties[0, 0].pressure.pprint()
    m.fs.unit.concentrate.properties[0, 1].pressure.pprint()
    m.fs.unit.friction_factor.pprint() if hasattr(m.fs.unit, "friction_factor") else 0
    print(
        value(m.fs.unit.diluate.properties[0, 0].pressure)
        - value(m.fs.unit.diluate.properties[0, 1].pressure)
    )
    m.fs.unit.hydraulic_diameter.pprint() if hasattr(
        m.fs.unit, "hydraulic_diameter"
    ) else 0
    m.fs.unit.visc_d.pprint()
    m.fs.unit.dens_mass.pprint()

    # m.fs.unit.diluate.properties[...].conc_mol_phase_comp["Liq","Na_+"].pprint()
    # m.fs.unit.diluate.properties[...].equiv_conductivity_phase["Liq"].pprint()
    # m.fs.unit.concentrate.properties[...].equiv_conductivity_phase["Liq"].pprint()

    if m.fs.unit.config.has_nonohmic_potential_membrane:
        m.fs.unit.conc_mem_surf_mol_x.display()
        m.fs.unit.potential_nonohm_membrane_x.display()
    if hasattr(m.fs.unit, "dl_thickness_x"):
        m.fs.unit.dl_thickness_x.pprint()
    if hasattr(m.fs.unit, "potential_ohm_dl_x"):
        m.fs.unit.potential_ohm_dl_x.pprint()
    if hasattr(m.fs.unit, "potential_nonohm_dl_x"):
        m.fs.unit.potential_nonohm_dl_x.pprint()
    """
    if m.fs.unit.config.diffusion_layer_polarization:
        for i in m.fs.unit.conc_mem_surf_mol_x:
            if i[4]=='Na_+':
                tx=(i[2], i[3])
                if (i[0] == "cem" and i[1] == "cathode_left") or (i[0] == "aem" and i[1] == "anode_right"):
                    rc = log(m.fs.unit.conc_mem_surf_mol_x[i].value * m.fs.unit.concentrate.properties[tx].conc_mol_phase_comp['Liq',i[4]].value **-1)
                    ri = log(1 + m.fs.unit.current_density_x[tx].value * m.fs.unit.current_dens_lim_x[tx].value **-1)
                    print(f'Ratio by conc: {i}, {rc}; Ratio by i: {ri}')
                else:
                    rc = log(m.fs.unit.conc_mem_surf_mol_x[i].value * m.fs.unit.diluate.properties[tx].conc_mol_phase_comp['Liq',i[4]].value **-1)
                    ri = log(1 - m.fs.unit.current_density_x[tx].value * m.fs.unit.current_dens_lim_x[tx].value **-1)
                    print(f'Ratio by conc: {i}, {rc}; Ratio by i: {ri}')
    """
    return m


if __name__ == "__main__":
    m = main()
