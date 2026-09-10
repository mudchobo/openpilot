"""
Copyright (c) 2021-, Haibin Wen, sunnypilot, and a number of other contributors.

This file is part of sunnypilot and is licensed under the MIT License.
See the LICENSE.md file in the root directory for more details.
"""
import pytest

from opendbc.car import structs

from openpilot.common.params import Params
from openpilot.common.test import OpenpilotTestCase
from openpilot.sunnypilot.selfdrive.controls.lib.latcontrol_torque_ext_override import \
  LatControlTorqueExtOverride, PARAM_READ_FRAMES, friction_scale

BASE_FRICTION = 0.1
BASE_LAT_ACCEL_FACTOR = 2.32


def _make_torque_params():
  CP = structs.CarParams()
  CP.lateralTuning.init('torque')
  CP.lateralTuning.torque.friction = BASE_FRICTION
  CP.lateralTuning.torque.latAccelFactor = BASE_LAT_ACCEL_FACTOR
  return CP.lateralTuning.torque


def _make_override(reduction=0, manual_override=False):
  params = Params()
  params.put_bool("EnforceTorqueControl", True, block=True)
  params.put("FrictionReduction", reduction, block=True)
  params.put_bool("TorqueParamsOverrideEnabled", manual_override, block=True)
  return LatControlTorqueExtOverride(structs.CarParams()), params


class TestFrictionScale(OpenpilotTestCase):
  def test_steps(self):
    assert friction_scale(0) == 1.0
    assert friction_scale(3) == pytest.approx(0.7)
    assert friction_scale(9) == pytest.approx(0.1)

  def test_clamped(self):
    assert friction_scale(-5) == 1.0
    assert friction_scale(100) == friction_scale(9)


class TestFrictionReduction(OpenpilotTestCase):
  def test_off_leaves_friction_untouched(self):
    ext, _ = _make_override(reduction=0)
    tp = _make_torque_params()
    for _ in range(PARAM_READ_FRAMES * 2):
      ext.update_override_torque_params(tp)
    assert tp.friction == pytest.approx(BASE_FRICTION)

  def test_does_not_compound_over_frames(self):
    """The scale is re-derived from the base each frame, never multiplied onto itself."""
    ext, _ = _make_override(reduction=3)
    tp = _make_torque_params()
    for _ in range(PARAM_READ_FRAMES * 2):
      ext.update_override_torque_params(tp)
    assert tp.friction == pytest.approx(BASE_FRICTION * 0.7)

  def test_follows_a_new_base_written_by_controlsd(self):
    """Self-Tune rewrites friction every frame; the scale must apply to the new value."""
    ext, _ = _make_override(reduction=2)
    tp = _make_torque_params()
    for _ in range(10):
      ext.update_override_torque_params(tp)
    assert tp.friction == pytest.approx(BASE_FRICTION * 0.8)

    learned = 0.15
    for _ in range(10):
      tp.friction = learned  # what controlsd does with the learned value
      ext.update_override_torque_params(tp)
    assert tp.friction == pytest.approx(learned * 0.8)

  def test_restores_base_when_turned_off(self):
    """With Self-Tune off nothing rewrites the base, so disarming has to restore it."""
    ext, params = _make_override(reduction=4)
    tp = _make_torque_params()
    for _ in range(PARAM_READ_FRAMES):
      ext.update_override_torque_params(tp)
    assert tp.friction == pytest.approx(BASE_FRICTION * 0.6)

    params.put("FrictionReduction", 0, block=True)
    for _ in range(PARAM_READ_FRAMES * 2):
      ext.update_override_torque_params(tp)
    assert tp.friction == pytest.approx(BASE_FRICTION)

  def test_lat_accel_factor_untouched(self):
    ext, _ = _make_override(reduction=5)
    tp = _make_torque_params()
    for _ in range(PARAM_READ_FRAMES * 2):
      ext.update_override_torque_params(tp)
    assert tp.latAccelFactor == pytest.approx(BASE_LAT_ACCEL_FACTOR)

  def test_no_limits_recalculation_requested(self):
    """friction doesn't feed the PID limits, so we never ask for update_limits()."""
    ext, _ = _make_override(reduction=3)
    tp = _make_torque_params()
    assert not any(ext.update_override_torque_params(tp) for _ in range(PARAM_READ_FRAMES * 2))

  def test_skipped_while_manual_override_active(self):
    ext, params = _make_override(reduction=5, manual_override=True)
    params.put("TorqueParamsOverrideFriction", 0.2, block=True)
    tp = _make_torque_params()
    for _ in range(PARAM_READ_FRAMES * 2):
      ext.update_override_torque_params(tp)
    assert tp.friction == pytest.approx(0.2)

  def test_disabled_without_enforce_torque_control(self):
    params = Params()
    params.put_bool("EnforceTorqueControl", False, block=True)
    params.put("FrictionReduction", 5, block=True)
    ext = LatControlTorqueExtOverride(structs.CarParams())
    tp = _make_torque_params()
    for _ in range(PARAM_READ_FRAMES * 2):
      ext.update_override_torque_params(tp)
    assert tp.friction == pytest.approx(BASE_FRICTION)
