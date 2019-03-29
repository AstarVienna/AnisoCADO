PSF Library
===========

A library of PSFs has been generated for use with SimCADO.

.. execute_code::
    :hide_code:

    import requests

    url = "https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/"
    file = "summary.txt"

    response = requests.get(url + file)
    response = response.text.split("\r\n")

    for line in response:
        print(line.split(".fits")[0])
        print("`{} <{}{}>.fits`_".format(line, url, line.split(".fits")[0]))

