import subprocess

try:
    import pyomo
except ImportError:
    print("Installing Pyomo")
    print(subprocess.run("pip install pyomo", shell=True))
    print(subprocess.run("conda install ipopt_bin -c cachemeorg", shell=True))
    print(subprocess.run("conda install glpk -c cachemeorg", shell=True))

try:
    import holoviews
except ImportError:
    print("Installing HoloViews")
    print(subprocess.run("conda install -c ioam holoviews bokeh", shell=True))
