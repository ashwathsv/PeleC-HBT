#include "prob.H"

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  {
    amrex::ParmParse pp("prob");
        pp.get("rho_0", PeleC::h_prob_parm_device->rho_0);
        pp.get("rho_r", PeleC::h_prob_parm_device->rho_r);
        pp.get("p_0", PeleC::h_prob_parm_device->p_0);
        pp.get("p_r", PeleC::h_prob_parm_device->p_r);
        pp.query("u_0", PeleC::h_prob_parm_device->u_0);
        pp.query("u_r", PeleC::h_prob_parm_device->u_r);
        pp.query("v_0", PeleC::h_prob_parm_device->v_0);
        pp.get("width", PeleC::h_prob_parm_device->width);
        pp.get("Mach_shock", PeleC::h_prob_parm_device->Mach_shock);
        pp.get("Mach_in", PeleC::h_prob_parm_device->Mach_in);
        pp.get("inflow_time", PeleC::h_prob_parm_device->inflow_time);
        pp.get("out_width", PeleC::h_prob_parm_device->out_width);
        pp.get("out_loend", PeleC::h_prob_parm_device->out_loend);
        pp.get("iprobe", PeleC::h_prob_parm_device->iprobe);
        pp.get("jprobe", PeleC::h_prob_parm_device->jprobe);
        pp.get("in_loend", PeleC::h_prob_parm_device->in_loend);

        pp.get("pressure_file", PeleC::h_prob_parm_device->pres_file);
  }

  amrex::Real L = (probhi[0] - problo[0]) * 0.2;

  amrex::Real cp = 0.0;
  PeleC::h_prob_parm_device->massfrac[0] = 1.0;
  auto eos = pele::physics::PhysicsType::eos();

  pele::physics::transport::TransParm trans_parm;

  // Default
  trans_parm.const_viscosity = 0.0;
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_conductivity = 0.0;
  trans_parm.const_diffusivity = 0.0;

#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(
    pele::physics::transport::trans_parm_g, &trans_parm, sizeof(trans_parm));
#else
  std::memcpy(
    pele::physics::transport::trans_parm_g, &trans_parm, sizeof(trans_parm));
#endif
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
