import pytest
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from electrodialysis_0D_Xiangyu import (
    ElectricalOperationMode,
    Electrodialysis0D,
)
from watertap.costing import WaterTAPCosting
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    value,
    Set,
    Param,
    Var,
    Constraint,
)
from idaes.core import (
    FlowsheetBlock,
    EnergyBalanceType,
    MaterialBalanceType,
    MomentumBalanceType,
)
from idaes.core import UnitModelCostingBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.util.testing import initialization_tester
