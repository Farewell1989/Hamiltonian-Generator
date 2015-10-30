from BasicClass.BaseSpacePy import *
def test_basespace():
    a=KSpace(reciprocals=[array([2*pi,0.0]),array([0.0,2*pi])],nk=100)
    a.plot(show=True)
    print a.volume/(2*pi)**2
    a=KSpace(reciprocals=[array([1.0,0.0]),array([0.5,sqrt(3.0)/2])],nk=100)
    a.plot(show=True)
    print a.volume
    square_gxm(nk=100).plot()

def test_basespace_functions():
    a=square_bz(reciprocals=[array([1.0,1.0]),array([1.0,-1.0])],nk=100)
    a.plot(show=True)
    print a.volume
    a=rectangle_bz(nk=100)
    a.plot(show=True)
    print a.volume/(2*pi)**2
    a=hexagon_bz(nk=100,vh='v')
    a.plot(show=True)
    print a.volume/(2*pi)**2
