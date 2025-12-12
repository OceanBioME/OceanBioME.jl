import OceanBioME.Particles: compute_buffer_variable

compute_buffer_variable(::Val{:Photosynthesis}, kelp::SugarKelp, u, v, w, T, NO₃, NH₄, PAR) =
    photosynthesis(kelp, T, PAR)
