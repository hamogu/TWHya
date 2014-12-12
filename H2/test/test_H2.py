from .. import H2

def test_H2numbers():
    Jup = 0
    Jlow = 1
    vup = 0
    vlow = 4
    
    assert H2.numbers2string(Jup, Jlow, vup, vlow) == '0-4 P(1)'
    assert H2.string2numbers('0-4 P(1)') == (Jup, Jlow, vup, vlow)
