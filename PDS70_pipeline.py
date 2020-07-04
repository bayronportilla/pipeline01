import cflux_alma
import cflux_jband
import profiles_modeled
import jband_ratios
import prepare_images
import azprofile_alma
import ratio_alma

# cflux_alma
image_alma=cflux_alma.image("RTout0001_000854.89.fits.gz",0.074,0.057,63.0)
cflux_alma.radial_profile(image_alma,120)


# cflux_jband
image_Qphi=cflux_jband.image("RToutObs0001_000001.25.fits.gz")
cflux_jband.radial_profile(image_Qphi,120.0)


# profiles_modeled
profiles_modeled.output_data()


# jband_ratios
jband_ratios.find_ratios()


# Prepare images
image_alma_rotated=prepare_images.prepare_alma_image(image_alma,158.6)
image_Qphi_rotated=prepare_images.prepare_Qphi_image(image_alma,158.6)
print(prepare_images.peak_flux_alma_model("alma_model_rotated.fits"))
print(prepare_images.peak_flux_Qphi_model("Qphi_model_rotated.fits"))


# azprofile
azprofile_alma.azimuthal_profile("alma_model_rotated.fits",0.004,70.0,30.0,49.7,158.6,113.43,72)


# Ratio alma
ratio_alma.ratio_SMsm()
