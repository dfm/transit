# -*- coding: utf-8 -*-

from __future__ import division, print_function

__all__ = ["test_impact", "test_period"]

import numpy as np
from ..transit import Body, Central, System


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
    try:
        body.b = 1.4
    except ValueError:
        pass
    else:
        assert False, "This impact parameter is too large. Test should fail"

    b0 = 0.2
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


def _measure_duration(body, delta=1e-6):
    body.t0 = 0.0

    # Compute the positive duration.
    t = np.arange(0, body.duration, delta)
    lc = body.system.light_curve(t)
    plus = t[1 - lc == 0][0]

    # Compute the negative duration.
    t = np.arange(0, -body.duration, -delta)
    lc = body.system.light_curve(t)
    minus = t[1 - lc == 0][0]

    return plus - minus


def test_duration():
    s = System(Central())
    body = Body(period=10, r=0.01, e=0.1)
    s.add_body(body)

    delta = np.abs(_measure_duration(body) - body.duration)
    print(delta)
    assert 0
