# -*- coding: utf-8 -*-

from __future__ import division, print_function

__all__ = ["test_impact", "test_period"]

import numpy as np
from .transit import Body, Central, System


def test_impact():
    s = System(Central())
    body = Body(period=10)
    s.add_body(body)

    # Basic conversion test.
    b0 = 0.5
    body.b = b0
    assert np.allclose(body.b, b0), "Basic impact parameter conversion failed"

    # Large impact parameter.
    b0 = 1.4
    body.b = b0
    assert np.allclose(body.b, b0), "Large impact parameter conversion failed"

    try:
        body.b = -0.5
    except ValueError:
        pass
    else:
        assert False, ("Negative impact parameters are invalid. "
                       "Test should fail.")

    # Eccentric orbit.
    body.e = 0.99
    b0 = 0.1
    body.b = b0
    assert np.allclose(body.b, b0), "Eccentric impact parameter failed"

    body.pomega = 0.5
    body.b = b0
    assert np.allclose(body.b, b0), ("Eccentric (with pomega) impact "
                                     "parameter failed")

    # What about with a mean observer angle?
    s = System(Central(), iobs=85.0)
    body = Body(period=10)
    s.add_body(body)
    b0 = 0.5
    body.b = b0
    assert np.allclose(body.b, b0), "Observer angle conversion failed"


def test_period():
    s = System(Central())
    body = Body(period=10)
    s.add_body(body)

    # Basic examples.
    for p0 in [0.1, 12.5, 5000.0, 1e5]:
        body.period = p0
        assert np.allclose(body.period, p0), \
            "Basic conversion failed for period {0}".format(p0)

    # Failure modes.
    try:
        body.period = -0.5
        body.period = 0.0
    except ValueError:
        pass
    else:
        assert False, ("Negative (or zero) periods are invalid. "
                       "Test should fail.")

    # Massive central body.
    for m0 in [0.1, 10.0, 5000.0]:
        s.central.mass = m0
        for p0 in [0.1, 12.5, 5000.0, 1e5]:
            body.period = p0
            assert np.allclose(body.period, p0), \
                "Massive central conversion failed for period {0}".format(p0)

    # Massive orbiting body.
    s.central.mass = 1.0
    for m0 in [0.1, 10.0, 5000.0]:
        body.mass = m0
        for p0 in [0.1, 12.5, 5000.0, 1e5]:
            body.period = p0
            assert np.allclose(body.period, p0), \
                "Massive body conversion failed for period {0}".format(p0)


def _measure_duration(nm, body, delta=5e-5):
    # Compute the positive duration.
    t = np.arange(body.t0, body.t0 + body.duration, delta)
    lc = body.system.light_curve(t)
    plus = t[1 - lc == 0][0]

    # Compute the negative duration.
    t = np.arange(body.t0, body.t0 - body.duration, -delta)
    lc = body.system.light_curve(t)
    minus = t[1 - lc == 0][0]

    print(plus - minus, body.duration)
    err = np.abs((plus - minus) - body.duration)
    assert err < 10*delta, "Duration test '{0}' failed.".format(nm)


def test_duration():
    s = System(Central())
    body = Body(period=10, r=0.01, b=0.5)
    s.add_body(body)

    # Basic tests.
    for p in [1.0, 10.0, 100.0, 1e4]:
        body.period = p
        body.b = 0.5
        _measure_duration("period = {0}".format(p), body)
    body.period = 10.0
