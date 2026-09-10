"""
Copyright (c) 2021-, Haibin Wen, sunnypilot, and a number of other contributors.

This file is part of sunnypilot and is licensed under the MIT License.
See the LICENSE.md file in the root directory for more details.
"""

from openpilot.common.params import Params

PARAM_READ_FRAMES = 300  # ~3s at 100Hz
MAX_FRICTION_REDUCTION = 9


def friction_scale(reduction: int) -> float:
  """Friction Reduction setting: step N lowers friction by N*10% (0 = off, 9 = -90%)."""
  return 1.0 - min(max(reduction, 0), MAX_FRICTION_REDUCTION) / 10.0


class LatControlTorqueExtOverride:
  def __init__(self, CP):
    self.CP = CP
    self.params = Params()
    self.enforce_torque_control_toggle = self.params.get_bool("EnforceTorqueControl")  # only during init
    self.torque_override_enabled = self.params.get_bool("TorqueParamsOverrideEnabled")
    self.friction_reduction = int(self.params.get("FrictionReduction", return_default=True))
    self.frame = -1

    # unscaled friction, and the scaled value we last wrote, so the scale is re-derived
    # instead of compounding across frames
    self._friction_base: float | None = None
    self._friction_written: float | None = None

  def _disarm_friction_reduction(self, torque_params) -> None:
    # nothing rewrites friction while Self-Tune is off, so put our scaling back before
    # letting go of it, otherwise the reduction sticks until reboot
    if self._friction_base is not None and torque_params.friction == self._friction_written:
      torque_params.friction = self._friction_base
    self._friction_base = None
    self._friction_written = None

  def update_friction_reduction(self, torque_params) -> None:
    """Scale the friction the controller uses down by the Friction Reduction setting.

    Runs every frame, right before get_friction() and torque_from_lateral_accel(). The base
    is whatever wrote torque_params.friction last: the offline value at init, or the learned
    value controlsd writes each frame while Self-Tune is on.
    """
    if self.friction_reduction == 0:
      self._disarm_friction_reduction(torque_params)
      return

    if torque_params.friction != self._friction_written:
      self._friction_base = torque_params.friction

    torque_params.friction = self._friction_base * friction_scale(self.friction_reduction)
    # read back rather than storing what we computed: torque_params.friction is a Float32,
    # so the stored value is rounded and a float64 copy would never compare equal
    self._friction_written = torque_params.friction

  def update_override_torque_params(self, torque_params) -> bool:
    if not self.enforce_torque_control_toggle:
      return False

    self.frame += 1
    read_params = self.frame % PARAM_READ_FRAMES == 0

    if read_params:
      self.torque_override_enabled = self.params.get_bool("TorqueParamsOverrideEnabled")
      self.friction_reduction = int(self.params.get("FrictionReduction", return_default=True))

    # Manual Real-Time Tuning sets friction to an exact value, don't scale it on top
    if self.torque_override_enabled:
      if read_params:
        torque_params.latAccelFactor = float(self.params.get("TorqueParamsOverrideLatAccelFactor", return_default=True))
        torque_params.friction = float(self.params.get("TorqueParamsOverrideFriction", return_default=True))
        self._friction_base = None
        self._friction_written = None
        return True
      return False

    self.update_friction_reduction(torque_params)
    return False  # friction doesn't feed the PID limits, no update_limits() needed
