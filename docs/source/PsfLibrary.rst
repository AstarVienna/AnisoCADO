PSF Library
===========

A library of PSFs has been generated for use with SimCADO.

..
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


`AnisoCADO_SCAO_FVPSF_4mas_EsoMedian_20190328.fits[2](0)-0.87um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/AnisoCADO_SCAO_FVPSF_4mas_EsoMedian_20190328.fits[2](0)-0.87um.png>`_

`AnisoCADO_SCAO_FVPSF_4mas_EsoMedian_20190328.fits[3](0)-1.05um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/AnisoCADO_SCAO_FVPSF_4mas_EsoMedian_20190328.fits[3](0)-1.05um.png>`_

`AnisoCADO_SCAO_FVPSF_4mas_EsoMedian_20190328.fits[4](0)-1.25um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/AnisoCADO_SCAO_FVPSF_4mas_EsoMedian_20190328.fits[4](0)-1.25um.png>`_

`AnisoCADO_SCAO_FVPSF_4mas_EsoMedian_20190328.fits[5](0)-1.65um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/AnisoCADO_SCAO_FVPSF_4mas_EsoMedian_20190328.fits[5](0)-1.65um.png>`_

`AnisoCADO_SCAO_FVPSF_4mas_EsoMedian_20190328.fits[6](0)-2.15um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/AnisoCADO_SCAO_FVPSF_4mas_EsoMedian_20190328.fits[6](0)-2.15um.png>`_

`AnisoCADO_SCAO_FVPSF_4mas_EsoQ1_20190328.fits[2](0)-0.87um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/AnisoCADO_SCAO_FVPSF_4mas_EsoQ1_20190328.fits[2](0)-0.87um.png>`_

`AnisoCADO_SCAO_FVPSF_4mas_EsoQ1_20190328.fits[3](0)-1.05um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/AnisoCADO_SCAO_FVPSF_4mas_EsoQ1_20190328.fits[3](0)-1.05um.png>`_

`AnisoCADO_SCAO_FVPSF_4mas_EsoQ1_20190328.fits[4](0)-1.25um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/AnisoCADO_SCAO_FVPSF_4mas_EsoQ1_20190328.fits[4](0)-1.25um.png>`_

`AnisoCADO_SCAO_FVPSF_4mas_EsoQ1_20190328.fits[5](0)-1.65um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/AnisoCADO_SCAO_FVPSF_4mas_EsoQ1_20190328.fits[5](0)-1.65um.png>`_

`AnisoCADO_SCAO_FVPSF_4mas_EsoQ1_20190328.fits[6](0)-2.15um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/AnisoCADO_SCAO_FVPSF_4mas_EsoQ1_20190328.fits[6](0)-2.15um.png>`_

`AnisoCADO_SCAO_FVPSF_4mas_EsoQ4_20190328.fits[2](0)-0.87um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/AnisoCADO_SCAO_FVPSF_4mas_EsoQ4_20190328.fits[2](0)-0.87um.png>`_

`AnisoCADO_SCAO_FVPSF_4mas_EsoQ4_20190328.fits[3](0)-1.05um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/AnisoCADO_SCAO_FVPSF_4mas_EsoQ4_20190328.fits[3](0)-1.05um.png>`_

`AnisoCADO_SCAO_FVPSF_4mas_EsoQ4_20190328.fits[4](0)-1.25um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/AnisoCADO_SCAO_FVPSF_4mas_EsoQ4_20190328.fits[4](0)-1.25um.png>`_

`AnisoCADO_SCAO_FVPSF_4mas_EsoQ4_20190328.fits[5](0)-1.65um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/AnisoCADO_SCAO_FVPSF_4mas_EsoQ4_20190328.fits[5](0)-1.65um.png>`_

`AnisoCADO_SCAO_FVPSF_4mas_EsoQ4_20190328.fits[6](0)-2.15um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/AnisoCADO_SCAO_FVPSF_4mas_EsoQ4_20190328.fits[6](0)-2.15um.png>`_

`Default_PSF_SCAO.fits[0](0)-1.65um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/Default_PSF_SCAO.fits[0](0)-1.65um.png>`_

`Default_PSF_SCAO.fits[1](0)-0.9um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/Default_PSF_SCAO.fits[1](0)-0.9um.png>`_

`Default_PSF_SCAO.fits[2](0)-1.2um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/Default_PSF_SCAO.fits[2](0)-1.2um.png>`_

`Default_PSF_SCAO.fits[3](0)-2.2um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/Default_PSF_SCAO.fits[3](0)-2.2um.png>`_

`MAORY_MCAO_FVPSF_1.5mas_20181203.fits[2](0)-1.635um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_MCAO_FVPSF_1.5mas_20181203.fits[2](0)-1.635um.png>`_

`MAORY_MCAO_FVPSF_1.5mas_20181203.fits[3](0)-0.86um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_MCAO_FVPSF_1.5mas_20181203.fits[3](0)-0.86um.png>`_

`MAORY_MCAO_FVPSF_1.5mas_20181203.fits[4](0)-1.245um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_MCAO_FVPSF_1.5mas_20181203.fits[4](0)-1.245um.png>`_

`MAORY_MCAO_FVPSF_1.5mas_20181203.fits[5](0)-2.145um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_MCAO_FVPSF_1.5mas_20181203.fits[5](0)-2.145um.png>`_

`MAORY_MCAO_FVPSF_4mas_20181203.fits[2](0)-1.635um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_MCAO_FVPSF_4mas_20181203.fits[2](0)-1.635um.png>`_

`MAORY_MCAO_FVPSF_4mas_20181203.fits[3](0)-0.86um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_MCAO_FVPSF_4mas_20181203.fits[3](0)-0.86um.png>`_

`MAORY_MCAO_FVPSF_4mas_20181203.fits[4](0)-1.245um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_MCAO_FVPSF_4mas_20181203.fits[4](0)-1.245um.png>`_

`MAORY_MCAO_FVPSF_4mas_20181203.fits[5](0)-2.145um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_MCAO_FVPSF_4mas_20181203.fits[5](0)-2.145um.png>`_

`MAORY_SCAO_FVPSF_1.5mas_20181203.fits[2](4)-1.635um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_SCAO_FVPSF_1.5mas_20181203.fits[2](4)-1.635um.png>`_

`MAORY_SCAO_FVPSF_1.5mas_20181203.fits[3](4)-0.86um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_SCAO_FVPSF_1.5mas_20181203.fits[3](4)-0.86um.png>`_

`MAORY_SCAO_FVPSF_1.5mas_20181203.fits[4](4)-1.245um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_SCAO_FVPSF_1.5mas_20181203.fits[4](4)-1.245um.png>`_

`MAORY_SCAO_FVPSF_1.5mas_20181203.fits[5](4)-2.145um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_SCAO_FVPSF_1.5mas_20181203.fits[5](4)-2.145um.png>`_

`MAORY_SCAO_FVPSF_1.5mas_20181203.fits[6](4)-1.02um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_SCAO_FVPSF_1.5mas_20181203.fits[6](4)-1.02um.png>`_

`MAORY_SCAO_FVPSF_4mas_20181203.fits[2](24)-1.635um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_SCAO_FVPSF_4mas_20181203.fits[2](24)-1.635um.png>`_

`MAORY_SCAO_FVPSF_4mas_20181203.fits[3](24)-0.86um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_SCAO_FVPSF_4mas_20181203.fits[3](24)-0.86um.png>`_

`MAORY_SCAO_FVPSF_4mas_20181203.fits[4](24)-1.245um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_SCAO_FVPSF_4mas_20181203.fits[4](24)-1.245um.png>`_

`MAORY_SCAO_FVPSF_4mas_20181203.fits[5](24)-2.145um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_SCAO_FVPSF_4mas_20181203.fits[5](24)-2.145um.png>`_

`MAORY_SCAO_FVPSF_4mas_20181203.fits[6](24)-1.02um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/MAORY_SCAO_FVPSF_4mas_20181203.fits[6](24)-1.02um.png>`_

`POPPY_NIR_Diff_limited_CONST_PSF_1mas_20181210.fits[0](0)-1.2um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/POPPY_NIR_Diff_limited_CONST_PSF_1mas_20181210.fits[0](0)-1.2um.png>`_

`POPPY_NIR_Diff_limited_CONST_PSF_1mas_20181210.fits[1](0)-1.6um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/POPPY_NIR_Diff_limited_CONST_PSF_1mas_20181210.fits[1](0)-1.6um.png>`_

`POPPY_NIR_Diff_limited_CONST_PSF_1mas_20181210.fits[2](0)-2.2um.png <https://www.univie.ac.at/simcado/InstPkgSvr/psfs/psf_summary/POPPY_NIR_Diff_limited_CONST_PSF_1mas_20181210.fits[2](0)-2.2um.png>`_

