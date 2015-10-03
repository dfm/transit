# -*- coding: utf-8 -*-

from __future__ import division, print_function

__all__ = ["test_impact", "test_period"]

import numpy as np
from .transit import Body, Central, System
from .simple import SimpleSystem


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
    dur = body.duration
    print(dur)
    t = np.arange(-2 * dur, 2 * dur, delta)
    lc = body.system.light_curve(t)
    measured = t[lc < 1.0].max() - t[lc < 1.0].min()
    assert np.allclose(measured, dur, atol=0.01), nm


def test_duration():
    s = System(Central())
    body = Body(period=10, r=0.01)
    s.add_body(body)
    for p in [1.0, 10.0, 100.0, 1e3]:
        body.period = p
        body.b = 0.5
        _measure_duration("period = {0}".format(p), body)


def _test_gradient(s, eps=1.345e-7, **kwargs):
    t = np.linspace(-1.0, 1.0, 5000)
    g = s.get_gradient(t, **kwargs)
    assert np.any(np.abs(g) > 0.0)

    vector = s.get_vector()
    names = s.get_parameter_names()
    for i, v in enumerate(vector):
        vector[i] = v + eps
        s.set_vector(vector)
        p = s.get_value(t, **kwargs)

        vector[i] = v - eps
        s.set_vector(vector)
        m = s.get_value(t, **kwargs)

        vector[i] = v
        s.set_vector(vector)

        n = names[i]
        assert np.allclose(0.5 * (p - m) / eps, g[:, i]), \
            "{0} ({1})".format(n, i)


def test_kepler_gradient(**kwargs):
    s = System(Central())
    s.add_body(Body(period=3.0, mass=0.1, radius=0.1, b=0.5, e=0.1, omega=0.5))
    s.thaw_parameter("*")
    _test_gradient(s, **kwargs)


def test_simple_gradient(**kwargs):
    s = SimpleSystem(period=3.0, ror=0.1, impact=0.5, duration=0.1)
    s.thaw_all_parameters()
    _test_gradient(s, **kwargs)


def test_modeling_protocol():
    from george.modeling import supports_modeling_protocol
    s = System(Central())
    s.add_body(Body(period=3.0, mass=0.1, radius=0.1, b=0.5, e=0.1, omega=0.5))
    assert supports_modeling_protocol(s, verbose=True)
