#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "mechanism.h"

#include "PeleC.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemDerive.H"
#include "prob_parm.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm)
{
  // Geometry

  state(i, j, k, URHO) = prob_parm.rho_0;
  state(i, j, k, UMX) = prob_parm.rho_0 * prob_parm.u_0;
  state(i, j, k, UMY) = prob_parm.rho_0 * prob_parm.v_0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = prob_parm.p_0 / (prob_parm.gamma-1.0);
  state(i, j, k, UEDEN) = state(i,j,k,UEINT) 
                        + prob_parm.rho_0 * (0.5 * (prob_parm.u_0 * prob_parm.u_0 +
                             prob_parm.v_0 * prob_parm.v_0));
  state(i, j, k, UTEMP) = prob_parm.p_0/(prob_parm.rho_0 * prob_parm.Rsp);
  for (int n = 0; n < NUM_SPECIES; n++)
    state(i, j, k, UFS + n) = prob_parm.rho_0 * prob_parm.massfrac[n];
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real* /*s_int[NVAR]*/,
  amrex::Real s_ext[NVAR],
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm)
{
  s_ext[URHO] = prob_parm.rho_0;
  s_ext[UMX] = prob_parm.rho_0 * prob_parm.u_0;
  s_ext[UMY] = prob_parm.rho_0 * prob_parm.v_0;
  s_ext[UMZ] = 0.0;
  s_ext[UEINT] = prob_parm.p_0 / (prob_parm.gamma-1.0);
  s_ext[UEDEN] = s_ext[UEINT] 
               + prob_parm.rho_0 *
                 (0.5 * (prob_parm.u_0 * prob_parm.u_0 +
                                          prob_parm.v_0 * prob_parm.v_0));
  s_ext[UTEMP] = prob_parm.p_0/(prob_parm.rho_0 * prob_parm.Rsp);
  for (int n = 0; n < NUM_SPECIES; n++)
    s_ext[UFS + n] = prob_parm.rho_0 * prob_parm.massfrac[n];
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_userbc(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NUM_GROW][NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int ind,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
    using namespace amrex;
    using amrex::Real;

    if(idir == 0 && sgn == -1){
        //supersonic inflow BC for left boundary
        Real ycm = Real(0.5)*(prob_parm.width + geomdata.ProbLo(1) + geomdata.ProbHi(1));
        Real tubelo = prob_parm.in_loend;
        Real tubehi = tubelo + prob_parm.width;
        // Print() << "ycm = " << ycm << ", tubelo= " << tubelo << ", tubehi= " << tubehi << "\n";
        if(x[1] >= tubelo && x[1] <= tubehi){
            if(time <= prob_parm.inflow_time){
                // high-pressure supersonic inflow BC
                s_ext[URHO]  = prob_parm.rho_r * prob_parm.rho_0;
                Real pre     = prob_parm.p_r * prob_parm.p_0;
                s_ext[UMY]   = Real(0.0);
                s_ext[UMZ]   = Real(0.0);
                s_ext[UEINT] = pre / (prob_parm.gamma - Real(1.0));
                Real cs      = std::sqrt(prob_parm.gamma * pre / (prob_parm.rho_r * prob_parm.rho_0));
                s_ext[UMX]   = prob_parm.rho_r * prob_parm.rho_0 * prob_parm.Mach_in * cs;
                s_ext[UEDEN] = s_ext[UEINT] 
                             + Real(0.5)*( s_ext[UMX]*s_ext[UMX] + s_ext[UMY]*s_ext[UMY] 
                                      + s_ext[UMZ]*s_ext[UMZ] )/s_ext[URHO];
                s_ext[UTEMP]  = pre/(prob_parm.Rsp * s_ext[URHO]);
                for (int nc = 0; nc < NUM_SPECIES; nc++)
                  s_ext[UFS + nc] = prob_parm.rho_0 * prob_parm.massfrac[nc];
            } else{
                // simple first-order extrapolation
                for(int nc = URHO; nc < NVAR; ++nc){
                    s_ext[nc] = s_int[0][nc];
                }
            }
        }else{
            // slip wall boundary condition (symmetry BC for scalars, no flux BC for normal velocity)
              if(prob_parm.isout_lox == 1 && x[1] >= geomdata.ProbHi(1) - prob_parm.width){
                // simple first-order extrapolation
                for(int nc = URHO; nc < NVAR; ++nc){
                    s_ext[nc] = s_int[0][nc];
                }
              }else{
                s_ext[URHO]   = s_int[-(ind+1)][URHO];
                s_ext[UMX]    = -s_int[-(ind+1)][UMX];
                s_ext[UMY]    = s_int[-(ind+1)][UMY];
                s_ext[UMZ]    = s_int[-(ind+1)][UMZ];
                s_ext[UEDEN]  = s_int[-(ind+1)][UEDEN];
                s_ext[UEINT]  = s_int[-(ind+1)][UEINT];
                s_ext[UTEMP]  = s_int[-(ind+1)][UTEMP];
                for (int nc = 0; nc < NUM_SPECIES; nc++)
                    s_ext[UFS + nc] = s_int[-(ind+1)][UFS+nc];
              }          
        }        
    }else if(idir == 1 && sgn == -1){
        // Outflow BC for Case 3 (outflow at bottom boundary)
        Real tubelo = prob_parm.out_loend;
        Real tubehi = tubelo + prob_parm.out_width;
        if(x[0] >= tubelo && x[0] <= tubehi){
            const amrex::Box& domain = geomdata.Domain();
            const auto domhi = amrex::ubound(domain);

            // simple zero-order extrapolation
            for(int nc = URHO; nc < NVAR; ++nc){
                s_ext[nc] = s_int[0][nc];
            }
      
        }else{
          // slip wall boundary condition (symmetry BC for scalars, no flux BC for normal velocity)
            s_ext[URHO]   = s_int[-(ind+1)][URHO];
            s_ext[UMX]    = s_int[-(ind+1)][UMX];
            s_ext[UMY]    = -s_int[-(ind+1)][UMY];
            s_ext[UMZ]    = s_int[-(ind+1)][UMZ];
            s_ext[UEDEN]  = s_int[-(ind+1)][UEDEN];
            s_ext[UEINT]  = s_int[-(ind+1)][UEINT];
            s_ext[UTEMP]  = s_int[-(ind+1)][UTEMP];
            for (int nc = 0; nc < NUM_SPECIES; nc++)
              s_ext[UFS + nc] = s_int[-(ind+1)][UFS+nc];
        } 
    }else if(idir == 1 && sgn == 1){
        // Outflow BC for Case 4 (outflow at top boundary)
        Real tubelo = prob_parm.out_loend;
        Real tubehi = tubelo + prob_parm.out_width;
        if(x[0] >= tubelo && x[0] <= tubehi){
            // simple first-order extrapolation
            for(int nc = URHO; nc < NVAR; ++nc){
                s_ext[nc] = s_int[0][nc];
            }
        }else{
            // slip wall boundary condition (symmetry BC for scalars, no flux BC for normal velocity)
            const amrex::Box& domain = geomdata.Domain();
            const auto domhi = amrex::ubound(domain);

            s_ext[URHO]   = s_int[ind-domhi.y-1][URHO];
            s_ext[UMX]    = s_int[ind-domhi.y-1][UMX];
            s_ext[UMY]    = -s_int[ind-domhi.y-1][UMY];
            s_ext[UMZ]    = s_int[ind-domhi.y-1][UMZ];
            s_ext[UEDEN]  = s_int[ind-domhi.y-1][UEDEN];
            s_ext[UEINT]  = s_int[ind-domhi.y-1][UEINT];
            s_ext[UTEMP]  = s_int[ind-domhi.y-1][UTEMP];
            for (int nc = 0; nc < NUM_SPECIES; nc++)
              s_ext[UFS+nc]  = s_int[ind-domhi.y-1][UFS+nc];
        }         
    }

}

// This function writes the pressure data from the probe to file
AMREX_GPU_DEVICE
inline
void
pc_probspecific_func (amrex::MultiFab& State, amrex::GeometryData const& geomdata, ProbParmDevice& prob_parm,
                        int stepnum, amrex::Real time)
{
    using namespace amrex;
    using amrex::Real;

    const Box& domain = geomdata.Domain();
    const auto domhi = amrex::ubound(domain);
    const auto domlo = amrex::lbound(domain);
    int myproc = ParallelDescriptor::MyProc();
    int tagproc = 0;

    Real probe_val[UFA] = {0.0};
    Real pres_val = 0.0;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(State,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& sfab = State.array(mfi);
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        // amrex::Print(myproc) << "myproc = " << myproc << ", lo = " << lo << ", hi = " << hi << "\n";

        if((prob_parm.iprobe >= lo.x && prob_parm.iprobe <= hi.x)
        && (prob_parm.jprobe >= lo.y && prob_parm.jprobe <= hi.y )){
            tagproc = myproc;
            amrex::Print(myproc) << "myproc = " << myproc << ", lo = " << lo << ", hi = " << hi 
                                 << "tagproc" << tagproc << "\n";
            for(int nc = URHO; nc < UFA; ++nc){
                probe_val[nc] = sfab(prob_parm.iprobe,prob_parm.jprobe,lo.z,nc);
            }
            // convert to SI units
            pres_val = probe_val[URHO] * prob_parm.Rsp * probe_val[UTEMP];
            probe_val[URHO]   = probe_val[URHO]*1000.0;
            probe_val[UMX]    = probe_val[UMX]*10.0;
            probe_val[UMY]    = probe_val[UMY]*10.0;
            probe_val[UMZ]    = probe_val[UMZ]*10.0;
            probe_val[UEINT]  = Real(0.1)*probe_val[UEINT];
            probe_val[UEDEN]  = Real(0.1)*probe_val[UEDEN];
            pres_val          = Real(0.1)*pres_val;
        }
    }

    if(myproc == tagproc){
        
        std::ofstream ofs;

        if(stepnum == 0){
            ofs.open(prob_parm.pres_file, std::ofstream::out);
            if(!ofs)  amrex::Abort("Failed to create file to write pressure data..exiting!");
            amrex::Print(myproc,ofs) << "# time rho rhou rhov rhow rhoE rhoe temp pre" << "\n";
            // Print(myproc) << "tagproc = " << tagproc << "\n";
        }else{
            ofs.open(prob_parm.pres_file, std::ofstream::app);
        }

        if(probe_val[URHO] > 0.0){
          amrex::Print(myproc, ofs).SetPrecision(8) << std::left << std::setw(12) << time << "\t"
        << std::left << std::setw(12) << probe_val[URHO]  << "\t"  
        << std::left << std::setw(12) << probe_val[UMX]   << "\t"
        << std::left << std::setw(12) << probe_val[UMY]   << "\t"
        << std::left << std::setw(12) << probe_val[UMZ]   << "\t"
        << std::left << std::setw(12) << probe_val[UEDEN] << "\t"
        << std::left << std::setw(12) << probe_val[UEINT] << "\t"
        << std::left << std::setw(12) << probe_val[UTEMP] << "\t"
        << std::left << std::setw(12) << pres_val  << "\n";
        }

        ofs.close();    
    }

    ParallelDescriptor::Barrier();

}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
