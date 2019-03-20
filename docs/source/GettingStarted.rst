Getting Started
===============

Determine the Strehl ratio for different wavelengths as we move away from the
(on-axis) guide star::

    sr_list = []
    wave_list = [0.8, 0.9, 1.1, 1.6, 2.15]
    for wave in wave_list:
        # Generate a PSF object for wavelength `wave` in um
        psf = AnalyticalScaoPsf(N=512, wavelength=wave)

        # Shift the PSF to different distances from the NGS
        sr = np.zeros(x.shape)
        for i in range(0, 41, 2):
            psf.shift_off_axis(x, 0)
            sr[i, j] = psf.strehl_ratio
        sr_list += [sr]

To plot the results::

    for wave, sr in zip(wave_list, sr_list):
        plt.plot(x.flatten(), sr, label="{} um".format(wave))
        # plt.semilogy()
    plt.legend()
    plt.xlabel("Distance from NGS [arcsec]")
    plt.ylabel("Strehl Ratio")
    plt.show()

.. figure:: ./_static/off_axis_strehl_ratios.png
    :align: center
    :figwidth: 700px
