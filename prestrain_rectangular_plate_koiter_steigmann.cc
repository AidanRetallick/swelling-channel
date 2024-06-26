// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
// LIC//
// LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
// LIC//
// LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
#include <fenv.h>
// Generic routines
#include "generic.h"

// The equations
#include "c1_koiter_steigmann.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;


template<class ELEMENT>
class MyMesh : public virtual TriangleMesh<ELEMENT>
{
public:
  MyMesh<ELEMENT>(TriangleMeshParameters& triangle_mesh_parameters,
                  TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
    : TriangleMesh<ELEMENT>(triangle_mesh_parameters, time_stepper_pt)
  {
  }

  void update_polyline_representation_from_restart()
  {
    oomph_info << "Called MyMesh::update_polyline_representation_from_restart()"
               << std::endl
               << "           NOT DOING ANYTHING IN THIS FUNCTION" << std::endl;
  }
};

namespace Parameters
{
  // Plate parameters
  double L1 = 10.0;
  double L2 = 1.0;
  double Initial_thickness = 0.00548;
  double Nu = 0.495;
  double Mu = 1.0e0;

  // Control parameters
  double P_mag = 0.0;
  double C_mag = 0.0;
  Data* C_swell_data_pt;
  double Isotropic_prestrain = 0.004; //0.00325;
  // Control parameter limits
  double P_max = 0.0;
  double C_swell_max = 0.0;
  double C_swell_min = 0.0;
  // Control parameter increment sizes
  double P_inc = 1.0;
  double C_swell_inc = 0.0;

  // Dimensions
  double L_dim = 4.93e-3;
  double E_dim = 2.03e6; //1.99e6;

  // Dependent parameters
  double Thickness = Initial_thickness * (1.0 + C_mag);
  double P_dim = E_dim * Initial_thickness;
  /// Function to update dependent parameters
  void update_dependent_parameters()
  {
    // Check the Pressure dim is correct
    P_dim = E_dim * Initial_thickness;
    // Update the swelling dependent data
    Parameters::C_swell_data_pt->set_value(0, Parameters::C_mag);
    Thickness = Initial_thickness * (1.0 + C_mag);
  }

  // Mesh parameters
  double element_area = 0.1;
  unsigned n_long_edge_nodes = 0;
  unsigned n_short_edge_nodes = 0;

  // Outputs
  string output_dir = "RESLT";
  ofstream* pvd_stream_pt;

  // Assigns the value of pressure depending on the position (x,y)
  void get_pressure(const Vector<double>& x, double& pressure)
  {
    // Constant pressure
    pressure = P_mag;
  }

  /// Pressure acting purely vertically
  void get_pressure_vert(const Vector<double>& x,
						 const Vector<double>& u,
			 const DenseMatrix<double>& grad_u,
			 const Vector<double>& n,
			 Vector<double>& pressure)
  {
    // Assign pressure
    pressure.resize(3);
    pressure[0] = 0.0;
    pressure[1] = 0.0;
    pressure[2] = P_mag;
  }

  /// Pressure depending on the position (x,y) and deformation of the sheet
  void get_pressure(const Vector<double>& x,
				const Vector<double>& u,
				const DenseMatrix<double>& grad_u,
				const Vector<double>& n,
				Vector<double>& pressure)
  {
    // Metric tensor of deformed surface
    DenseMatrix<double> G(2,2,0.0);
    for(unsigned alpha = 0; alpha < 2; alpha++)
    {
      G(alpha, alpha) += 1.0;
      for (unsigned beta = 0; beta < 2; beta++)
      {
	G(alpha, beta) += grad_u(alpha, beta) + grad_u(beta, alpha);
	for (unsigned i = 0; i < 3; i++)
	{
		G(alpha, beta) += grad_u(i, alpha) * grad_u(i, beta);
	}
      }
    }

    // Find the pressure per undeformed area in terms of the pressure per
    // deformed area
    double p = sqrt(G(0,0)*G(1,1) - G(1,0)*G(0,1)) * P_mag;

    // Assign pressure
    pressure.resize(3);
    pressure[0] = p * n[0];
    pressure[1] = p * n[1];
    pressure[2] = p * n[2];
  }

  // Assigns the value of swelling depending on the position (x,y)
  void get_swelling(const Vector<double>& x, double& swelling)
  {
    // Almost uniform swelling with clamped boundaries
    swelling =
      C_swell_data_pt->value(0); // * (1 - pow(x[0]/L1,10) - pow(x[1]/L2,10)
    // + pow( (x[0]*x[1])/(L1*L2) , 10));
  }


  /// Swelling induced prestrain
  void get_swelling_prestrain(const Vector<double>& x,
						DenseMatrix<double>& prestrain)
  {
    // Swelling at x
    double c = 0.0;
    get_swelling(x, c);
    double isostrain = Isotropic_prestrain - (c + 0.5*c*c);
    prestrain(0,0) = isostrain;
    prestrain(0,1) = 0.0;
    prestrain(1,0) = 0.0;
    prestrain(1,1) = isostrain;
  }


  //----------------------------------------------------------------------
  // Mooney-Rivlin stress
  //----------------------------------------------------------------------
  // According to Li & Healey (2016)
  //   C1 + C2 = E / 6.0
  // Therefore, non-dimensionalising by E, we get that
  //   C2 = 1.0 / 6.0 - C1

  /// First dimensionless Mooney-Rivlin constant
  //double C1 = 1.0 / 6.6;
  double C1 = 0.1994;

  // Calculate in stress to ensure both don't fall out of sync
  // /// Second dimensionless Mooney-Rivlin constant
  // double C2 = 1.0 / 6.0 - C1;

  /// Mooney Rivlin stress function
  void mooney_rivlin_stress(const Vector<double>& x,
			    const Vector<double>& u,
			    const DenseMatrix<double>& e,
			    const DenseMatrix<double>& g,
			    DenseMatrix<double>& stress)
  {
    // Constants
    const double c1 = Parameters::C1;
    const double c2 = 1.0 / 6.0 - c1;

    // Matrix of cofactors of strain tensor
    DenseMatrix<double> cof_e(2, 2);
    cof_e(0, 0) = e(1, 1);
    cof_e(1, 1) = e(0, 0);
    cof_e(0, 1) = -e(0, 1);
    cof_e(1, 0) = -e(1, 0);
    // Matrix of cofactors of metric tensor
    DenseMatrix<double> cof_g(2, 2);
    cof_g(0, 0) = g(1, 1);
    cof_g(1, 1) = g(0, 0);
    cof_g(0, 1) = -g(0, 1);
    cof_g(1, 0) = -g(1, 0);

    // Determinants
    const double det_e = e(0, 0) * e(1, 1) - e(0, 1) * e(1, 0);
    const double det_g = g(0, 0) * g(1, 1) - g(0, 1) * g(1, 0);
    // Traces
    const double tr_e = e(0, 0) + e(1, 1);
    // const double tr_g = g(0,0)+g(1,1);
    // NB det(g) = 4 det(e) + 2 Tr(e) +1
    // Determinant of g squared minus one
    const double det_g_2_m1 =
      (4 * det_e + 2 * tr_e) * (4 * det_e + 2 * tr_e + 2);

    // Now fill in the stress
    // Loop over indices
    DenseMatrix<double> i2(2, 2, 0.0);
    i2(0, 0) = 1.0;
    i2(1, 1) = 1.0;

    // Now Fill in the Stress
    for (unsigned alpha = 0; alpha < 2; ++alpha)
    {
      for (unsigned beta = 0; beta < 2; ++beta)
      {
        // 2nd Piola Kirchhoff (Membrane) Stress for Mooney Rivlin model
        //  stress(alpha,beta)=2*(c1+ c2 / det_g) * kronecker(alpha,beta)
        //       + 2*((- c1 - tr_g *c2) / pow(det_g,2) + c2)*cof_g(alpha,beta);
        // For 2D:
        // Cof g = I + 2 Cof e
        // tr(g) = 2 + 2 tr(e)
        // Det(g) = 4 Det(e) + 2 Tr(e) + 1
        // 2nd Piola Kirchhoff (Membrane) Stress for Mooney Rivlin model
        // ( Modified so that all terms are order epsilon if possible *)
        stress(alpha, beta) =
	  2 * (c2 * (-4 * det_e - 2 * tr_e) / det_g) * i2(alpha, beta) -
	  4 * ((c1 + 2 * c2 + 2 * tr_e * c2) / pow(det_g, 2) - c2) *
	    cof_e(alpha, beta) -
          2 *
            ((c1 + 2 * c2 + 2 * tr_e * c2) * (-det_g_2_m1) / pow(det_g, 2) +
             2 * tr_e * c2) *
            i2(alpha, beta);
      }
    }
  }

  /// Mooney Rivlin stiffness tensor (only fills in dstrain, not du)
  void d_mooney_rivlin_stress_d_strain(const Vector<double>& x,
				       const Vector<double>& u,
				       const DenseMatrix<double>& strain,
				       const DenseMatrix<double>& g,
                                       RankThreeTensor<double>& d_stress_du,
                                       RankFourTensor<double>& d_stress_dstrain)
  {
    // Constants
    const double c1 = Parameters::C1;
    const double c2 = 1.0 / 6.0 - c1;

    // Matrix of cofactors of metric tensor
    DenseMatrix<double> cof_g(2, 2);
    cof_g(0, 0) = g(1, 1);
    cof_g(1, 1) = g(0, 0);
    cof_g(0, 1) = -g(0, 1);
    cof_g(1, 0) = -g(1, 0);

    // Fill in determinants
    const double det_g = g(0, 0) * g(1, 1) - g(0, 1) * g(1, 0);
    const double tr_g = g(0, 0) + g(1, 1);

    // Identity matrix
    DenseMatrix<double> i2(2, 2, 0.0);
    i2(0, 0) = 1.0;
    i2(1, 1) = 1.0;

    // Now Fill in the Stress
    for (unsigned alpha = 0; alpha < 2; ++alpha)
    {
      for (unsigned beta = 0; beta < 2; ++beta)
      {
        // 2nd Piola Kirchhoff (Membrane) Stress for Mooney Rivlin model
        // stress(alpha,beta)=2*(c1+ c2 / det_g) * i2(alpha,beta)
        //      + 2*((- c1 - tr_g *c2) / pow(det_g,2) + c2)*cof_g(alpha,beta);
        for (unsigned gamma = 0; gamma < 2; ++gamma)
        {
          for (unsigned delta = 0; delta < 2; ++delta)
          {
            // 2nd Piola Kirchhoff (Membrane) Stress for Mooney Rivlin model
            // d (cof_g) dg
            d_stress_dstrain(alpha, beta, gamma, delta) =
              2 * (-2 * (c2 / pow(det_g, 2)) * i2(alpha, beta) *
                     cof_g(gamma, delta) +
                   2 * (c2 - (c1 + tr_g * c2) / pow(det_g, 2)) *
                     // dcof(g)/dg = (cof(g) \otimes cof(g) - cof(g) . dg / dg .
                     // cof(g))/det(g)
                     (cof_g(alpha, beta) * cof_g(gamma, delta) -
                      cof_g(alpha, gamma) * cof_g(delta, beta)) /
                     det_g +
                   4 * ((c1 + tr_g * c2) / pow(det_g, 3)) * cof_g(alpha, beta) *
                     cof_g(gamma, delta) -
                   2 * (c2 / pow(det_g, 2)) * cof_g(alpha, beta) *
                     i2(gamma, delta));
          }
        }
      }
    }
  }



  // Get the exact solution
  void get_null_fct(const Vector<double>& X, double& exact_w)
  {
    exact_w = 0.0;
  }

  // Get the exact solution(s)
  void dummy_exact_w(const Vector<double>& x, Vector<double>& exact_w)
  {
    /* Do nothing -> no exact solution in this case */
  }

} // namespace Parameters

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class UnstructuredKSProblem : public virtual Problem
{
public:
  /// Constructor
  UnstructuredKSProblem();

  /// Destructor
  ~UnstructuredKSProblem()
  {
    Trace_file_dim.close();
    Trace_file_nondim.close();
    delete (Surface_mesh_pt);
    delete (Bulk_mesh_pt);
    // Clean up memory
    delete Boundary_pt;
    //  delete Inner_open_boundaries_pt[0];
    //  delete Inner_open_boundaries_pt[1];
    delete Boundary0_pt;
    delete Boundary1_pt;
    delete Boundary2_pt;
    delete Boundary3_pt;
  };

  /// Setup and build the mesh
  void build_mesh();

  /// Return the number of newton iterations in the last solve attempt
  unsigned nnewton_iter_taken()
  {
    return Nnewton_iter_taken;
  }

  /// Temporal error norm function.
  double global_temporal_error_norm()
  {
    double global_error = 0.0;

    // Find out how many nodes there are in the problem
    unsigned n_node = mesh_pt()->nnode();

    // Loop over the nodes and calculate the estimated error in the values
    for (unsigned i = 0; i < n_node; i++)
    {
      double error = 0;
      Node* node_pt = mesh_pt()->node_pt(i);
      // Get error in solution: Difference between predicted and actual
      // value for nodal value 2 (only if we are at a vertex node)
      if (node_pt->nvalue() > 2)
      {
        error = node_pt->time_stepper_pt()->temporal_error_in_value(
          mesh_pt()->node_pt(i), 2);
      }
      // Add the square of the individual error to the global error
      global_error += error * error;
    }

    // Divide by the number of nodes
    global_error /= double(n_node);

    // Return square root...
    return sqrt(global_error);
  }

  void enable_xz_profiles()
  {
    STORE_XZ_PROFILES = true;
  }

  void disable_xz_profiles()
  {
    STORE_XZ_PROFILES = false;
  }

  /// Update after solve (empty)
  void actions_after_newton_solve()
  {
    /* No actions before newton solve */
  }

  /// Pin the in-plane displacements and set to zero at centre
  void pin_in_plane_displacements_at_centre_node();

  // /// Things to repeat after every newton iteration
  // void actions_before_newton_step()
  // {
  //   // File prefix and suffix strings
  //   std::string prefix = Doc_unsteady_info.directory() + "/";
  //   std::string suffix = std::to_string(Doc_unsteady_info.number()) + ".txt";

  //   // Get the residual and jacobian
  //   LinearAlgebraDistribution* dist = dof_distribution_pt();
  //   DoubleVector residual(dist);
  //   CRDoubleMatrix jacobian(dist);
  //   get_jacobian(residual, jacobian);
  //   residual.output(prefix + "residual" + suffix);
  //   jacobian.sparse_indexed_output(prefix + "jacobian" + suffix);
  // }

  /// Print information about the parameters we are trying to solve for.
  void actions_before_newton_solve()
  {
    // Print some gubbins about parameters
    oomph_info << "-------------------------------------------------------"
               << std::endl;
    oomph_info << "Solving for p = " << Parameters::P_mag << "  ("
               << Parameters::P_mag * Parameters::P_dim << "Pa)" << std::endl;
    oomph_info << "      c_swell = " << Parameters::C_swell_data_pt->value(0)
               << std::endl;
    oomph_info << "     Doc_info = " << Doc_steady_info.number() << ", "
               << Doc_unsteady_info.number() << std::endl;
    oomph_info << "-------------------------------------------------------"
               << std::endl;
  }

  /// Remove surface mesh before reading
  void actions_before_read_unstructured_meshes()
  {
    // How many surface elements are in the surface mesh
    unsigned n_element = Surface_mesh_pt->nelement();

    // Loop over the surface elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Kill surface element
      delete Surface_mesh_pt->element_pt(e);
    }

    Surface_mesh_pt->flush_element_and_node_storage();

    rebuild_global_mesh();
  }

  /// Add surface mesh back after reading
  void actions_after_read_unstructured_meshes()
  {
    rebuild_global_mesh();
    apply_boundary_conditions();
    complete_problem_setup();
  }


  /// Used damped solves to get close to a steady solution, when close
  /// enough, attempt a steady solve. If that fails, be stricter about the
  /// meaning of "close enough" and repeat until a steady solve succeeds.
  ///
  /// Expects:
  ///   dt_supplied_guess -- a guess for a good timestep size
  ///             epsilon -- an 'error tolerance' for the timestepper to limit
  ///                        the size of a damped step
  /// Accepts:
  ///        doc_unsteady -- do we document every unsteady (nonphysical) solve
  ///   begin_with_steady -- do we begin by trying a steady solve
  ///
  /// Returns:
  ///   A tuple containing suggestions to use as the inputs in the next
  ///   damped_solve. First a double suggestion for dt which is the timestep
  ///   of the first successful solve -- if the size of deformation is roughly
  ///   the same, it should work fine. Secondly a bool suggestion for
  ///   begin_with_steady which is true if this damped_solve took <2 timesteps.
  std::tuple<double, bool> damped_solve(const double& dt,
					const double& epsilon,
					const bool& doc_unsteady = false,
					const bool& begin_with_steady = false);

  /// Doc the solution
  void doc_solution(const bool steady, const std::string& comment = "");

  /// Dump problem to disk to allow for restart.
  void dump_it(ofstream& dump_file)
  {
    dump_file << Doc_steady_info.number() << " # steady step number"
              << std::endl
              << Doc_unsteady_info.number() << " # unsteady step number"
              << std::endl
              << Parameters::P_mag << " # pressure" << std::endl
              << Parameters::C_swell_data_pt->value(0) << " # swelling" << std::endl
              << Next_dt << " # suggested next timestep" << std::endl;

    // Call generic dump()
    Problem::dump(dump_file);
  }

  /// Read problem for restart from specified restart file.
  void restart(ifstream& restart_file)
  {
    string input_string;

    // Read line up to termination sign
    getline(restart_file, input_string, '#');
    // Ignore rest of line
    restart_file.ignore(80, '\n');
    // Read in steady step number
    Doc_steady_info.number() = unsigned(atof(input_string.c_str()));
    Doc_steady_info.number()++;

    getline(restart_file, input_string, '#');
    // Ignore rest of line
    restart_file.ignore(80, '\n');
    // Read in unsteady step number
    Doc_unsteady_info.number() = unsigned(atof(input_string.c_str()));
    Doc_unsteady_info.number() = 0; // for now [witnessme]

    getline(restart_file, input_string, '#');
    // Ignore rest of line
    restart_file.ignore(80, '\n');
    // Read in pressure value
    Parameters::P_mag = double(atof(input_string.c_str()));

    getline(restart_file, input_string, '#');
    // Ignore rest of line
    restart_file.ignore(80, '\n');
    // Read in steady step number
    Parameters::C_mag = double(atof(input_string.c_str()));
    Parameters::update_dependent_parameters();

    // Read line up to termination sign
    getline(restart_file, input_string, '#');
    // Ignore rest of line
    restart_file.ignore(80, '\n');
    // Read suggested next timestep
    Next_dt = double(atof(input_string.c_str()));

    // Refine the mesh and read in the generic problem data
    Problem::read(restart_file);

  } // end of restart

  /// \short Overloaded version of the problem's access function to
  /// the mesh. Recasts the pointer to the base Mesh object to
  /// the actual mesh type.
  MyMesh<ELEMENT>* mesh_pt()
  {
    return Bulk_mesh_pt;
  }

  /// Access function for Next_dt
  double next_dt()
  {
    return Next_dt;
  }

  /// Doc info object for labeling all output
  DocInfo Doc_steady_info;
  DocInfo Doc_unsteady_info;

private:
  // Triangle Mesh Parameter Data
  // This is the data used to set-up the mesh, we need to store the pointers
  // HERE otherwise we will not be able to clean up the memory once we have
  // finished the problem.
  // The initial (and maximum) element area
  double Element_area;
  // The mesh parameters
  TriangleMeshParameters* Triangle_mesh_parameters_pt;
  TriangleMeshClosedCurve* Boundary_pt;
  Vector<TriangleMeshOpenCurve*> Inner_open_boundaries_pt;
  TriangleMeshPolyLine* Boundary0_pt;
  TriangleMeshPolyLine* Boundary1_pt;
  TriangleMeshPolyLine* Boundary2_pt;
  TriangleMeshPolyLine* Boundary3_pt;
  TriangleMeshPolyLine* InnerBoudary0_pt;
  TriangleMeshPolyLine* InnerBoudary1_pt;

  /// Helper function to apply boundary conditions
  void apply_boundary_conditions();

  /// \short Helper function to (re-)set boundary condition
  /// and complete the build of all elements
  void complete_problem_setup();

  /// Trace file to document norm of solution
  ofstream Trace_file_dim, Trace_file_nondim;

  // Keep track of boundary ids
  enum
  {
    Outer_boundary0 = 0,
    Outer_boundary1 = 1,
    Outer_boundary2 = 2,
    Outer_boundary3 = 3
  };

  /// \short Delete traction elements and wipe the surface mesh
  void delete_traction_elements(Mesh* const& surface_mesh_pt);

  /// Pointer to "bulk" mesh
  MyMesh<ELEMENT>* Bulk_mesh_pt;

  /// Pointer to "surface" mesh
  Mesh* Surface_mesh_pt;

  /// Suggestion for next timestep (stored to allow it to be written
  /// to (or read from) restart file)
  double Next_dt = 1.0;

  /// Do we want to locate_zeta along the centreline and store the profile?
  bool STORE_XZ_PROFILES = false;

}; // end_of_problem_class

/// Constructor definition
template<class ELEMENT>
UnstructuredKSProblem<ELEMENT>::UnstructuredKSProblem()
  : Element_area(Parameters::element_area)
{
  add_time_stepper_pt(new BDF<2>(true));

  Problem::Always_take_one_newton_step = true;
  // Build the mesh
  build_mesh();

  // Store number of bulk elements
  complete_problem_setup();

  char filename[100];
  ofstream Param_file;
  strcpy(filename, (Parameters::output_dir + "/parameters.dat").c_str());
  Param_file.open(filename);

  // Output plate parameters
  Param_file << "L1           " << Parameters::L1 << std::endl
             << "L2           " << Parameters::L2 << std::endl
             << "thickness    " << Parameters::Thickness << std::endl
             << "nu           " << Parameters::Nu << std::endl
             << std::endl
             << "L_dim        " << Parameters::L_dim << std::endl
             << "E_dim        " << Parameters::E_dim << std::endl
             << "P_dim        " << Parameters::P_dim << std::endl
             << std::endl
             << "Element area " << Parameters::element_area << std::endl
             << "N_xedgenodes " << Parameters::n_long_edge_nodes << std::endl
             << "N_yedgenodes " << Parameters::n_short_edge_nodes << std::endl;

  strcpy(filename, (Parameters::output_dir + "/trace_nondim.dat").c_str());
  Trace_file_nondim.open(filename);
  strcpy(filename, (Parameters::output_dir + "/trace_dim.dat").c_str());
  Trace_file_dim.open(filename);

#ifdef OOMPH_HAS_MUMPS
  // Let everyone know we are going to use MUMPS
  oomph_info << std::endl << "Using MUMPS solver" << std::endl << std::endl;

  // Change solver
  linear_solver_pt() = new MumpsSolver;

  // Shut up
  dynamic_cast<MumpsSolver*>(linear_solver_pt())
    ->enable_suppress_warning_about_MPI_COMM_WORLD();

#endif

  oomph_info << "Number of equations: " << assign_eqn_numbers() << '\n';

  // Output the mesh
  strcpy(filename, (Parameters::output_dir + "/mesh.dat").c_str());
  ofstream meshfile(filename);
  Problem::mesh_pt()->output(meshfile);
  meshfile.close();

} // end Constructor

/// Set up and build the mesh
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::build_mesh()
{
  // Rectangular mesh boundary
  //
  //   V3        E2       V2
  //    O--------o--------O     ^
  //    |        V6       |     |
  //    |                 |     |
  // E3 oV7      oV8    V5o E1  L2
  //    |                 |     |
  //    |        V4       |     |
  //    O--------o--------O     v
  //   V0        E0       V1
  //    <--------L1------->

  double L1 = Parameters::L1;
  double L2 = Parameters::L2;


  unsigned Nl = Parameters::n_long_edge_nodes;
  unsigned Ns = Parameters::n_short_edge_nodes;
  double hl = L1 / (Nl - 1);
  double hs = L2 / (Ns - 1);
  Vector<Vector<double>> E0(Nl, Vector<double>(2, 0.0)),
    E1(Ns, Vector<double>(2, 0.0)), E2(Nl, Vector<double>(2, 0.0)),
    E3(Ns, Vector<double>(2, 0.0));
  for (unsigned i = 0; i < Nl; i++)
  {
    E0[i][0] = -0.5 * L1 + i * hl;
    E0[i][1] = -0.5 * L2;
    E2[i][0] = 0.5 * L1 - i * hl;
    E2[i][1] = 0.5 * L2;
  }
  for (unsigned i = 0; i < Ns; i++)
  {
    E1[i][0] = 0.5 * L1;
    E1[i][1] = -0.5 * L2 + i * hs;
    E3[i][0] = -0.5 * L1;
    E3[i][1] = 0.5 * L2 - i * hs;
  }

  // Define boundaries from edges.
  Boundary0_pt = new TriangleMeshPolyLine(E0, 0);
  Boundary1_pt = new TriangleMeshPolyLine(E1, 1);
  Boundary2_pt = new TriangleMeshPolyLine(E2, 2);
  Boundary3_pt = new TriangleMeshPolyLine(E3, 3);

  Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);
  boundary_polyline_pt[0] = Boundary0_pt;
  boundary_polyline_pt[1] = Boundary1_pt;
  boundary_polyline_pt[2] = Boundary2_pt;
  boundary_polyline_pt[3] = Boundary3_pt;

  Boundary_pt = new TriangleMeshClosedCurve(boundary_polyline_pt);

  TriangleMeshParameters Triangle_mesh_parameters(Boundary_pt);

  // Set the maximum element area
  Triangle_mesh_parameters.element_area() = Element_area;

  // Build an assign bulk mesh
  Bulk_mesh_pt =
    new MyMesh<ELEMENT>(Triangle_mesh_parameters, time_stepper_pt());

  // Create "surface mesh" that will contain only the prescribed-traction
  // elements. The constructor creates the mesh without adding any nodes
  // elements etc.
  Surface_mesh_pt = new Mesh;

  // Add two submeshes to problem
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Surface_mesh_pt);

  // Combine submeshes into a single Mesh
  build_global_mesh();
} // end build_mesh


//==start_of_complete======================================================
/// Set boundary condition exactly, and complete the build of
/// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::complete_problem_setup()
{
  DTSF_max_increase = 2.0;
  DTSF_min_decrease = 0.5;

  // Create a new Data object whose one-and-only value contains the
  // (in principle) adjustible load
  Parameters::C_swell_data_pt = new Data(1);

  // And pin the pressure as it is our control.
  Parameters::C_swell_data_pt->pin(0);
  Parameters::C_swell_data_pt->set_value(0, Parameters::C_mag);

  // Complete the build of all elements so they are fully functional
  unsigned n_element = Bulk_mesh_pt->nelement();
  for (unsigned e = 0; e < n_element; e++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    // Set the pressure & prestrain function pointers
    el_pt->pressure_fct_pt() = &Parameters::get_pressure;
    el_pt->prestrain_fct_pt() = &Parameters::get_swelling_prestrain;

    // Assign the parameter pointers for the element
    el_pt->thickness_pt() = &Parameters::Thickness;
    el_pt->nu_pt() = &Parameters::Nu;
    el_pt->mu_pt() = &Parameters::Mu;

    // Enable damping in all-direction
    el_pt->enable_damping();

    if(!CommandLineArgs::command_line_flag_has_been_set("--use_linear_stress"))
    {
      el_pt->stress_fct_pt() = &Parameters::mooney_rivlin_stress;
      el_pt->d_stress_fct_pt() = &Parameters::d_mooney_rivlin_stress_d_strain;
    }
    if(CommandLineArgs::command_line_flag_has_been_set("--use_fd_jacobian"))
    {
      el_pt->enable_finite_difference_jacobian();
    }

  }

  // Do we want to pin in-plane displacement?
  if (CommandLineArgs::command_line_flag_has_been_set("--pininplane"))
  {
    cout << "gonna pin em" << endl;
    // Pin the in-plane displacements
    unsigned nnode = mesh_pt()->nnode();
    for (unsigned inode = 0; inode < nnode; inode++)
    {
      mesh_pt()->node_pt(inode)->pin(0);
      mesh_pt()->node_pt(inode)->pin(1);
    }
    cout << "successfully pinned" << endl;
    assign_eqn_numbers();
  }

  // Set the boundary conditions
  apply_boundary_conditions();
}



//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::apply_boundary_conditions()
{
  // Displacement dofs:
  //-------------------
  // |   0   |   1   |   2   |   3   |   4   |   5   |
  // |  u_i  | u_i_x | u_i_y | u_i_xx| u_i_xy| u_i_yy|

  // Case: The plate is pinned (w given, dw/dn left free) along a boundary where
  // the outer unit normal points in the postive or negative x direction, so x
  // is constant and y varies along the boundary.  We therefore have to pin (and
  // assign values for) w, dw/dy and d^2w/dy^2
  const Vector<unsigned> pinned_edge_xn_dof{0, 2, 5};

  // Case: The plate is pinned (w given, dw/dn left free) along a boundary where
  // the outer unit normal points in the postive or negative y direction, so y
  // is constant and x varies along the boundary.  We therefore have to pin (and
  // assign values for) w, dw/dx and d^2w/dx^2
  const Vector<unsigned> pinned_edge_yn_dof{0, 1, 3};

  // Case: The plate is sliding (w left free, dw/dn given) along a boundary
  // where the outer unit normal points in the postive or negative x direction,
  // so x is constant and y varies along the boundary.  We therefore have to pin
  // (and assign values for) dw/dx and d^2w/dxdy
  const Vector<unsigned> sliding_clamp_xn_dof{1, 4};

  // Case: The plate is sliding (w left free, dw/dn given) along a boundary
  // where the outer unit normal points in the postive or negative y direction,
  // so y is constant and x varies along the boundary.  We therefore have to pin
  // (and assign values for) dw/dy and d^2w/dxdy
  const Vector<unsigned> sliding_clamp_yn_dof{2, 4};

  // Case: The plate is clamped (w given, dw/dn given) along a boundary where
  // the outer unit normal points in the postive or negative x direction, so x
  // is constant and y varies along the boundary.  We therefore have to pin (and
  // assign values for) w, dw/dx, dw/dy, d^2w/dxdy and d^2w/dy^2
  const Vector<unsigned> fully_clamped_xn_dof{0, 1, 2, 4, 5};

  // Case: The plate is clamped (w given, dw/dn given) along a boundary where
  // the outer unit normal points in the postive or negative y direction, so y
  // is constant and x varies along the boundary.  We therefore have to pin (and
  // assign values for) w, dw/dx, dw/dy, d^2w/dx^2 and d^2w/dxdy
  const Vector<unsigned> fully_clamped_yn_dof{0, 1, 2, 3, 4};

  // Free edge has no constraints
  const Vector<unsigned> free{};

  // [zdec] NONPHYSICAL - used for debugging
  const Vector<unsigned> uber_clamped{0, 1, 2, 3, 4, 5};

  //------------------------------------------------------------------
  //------------------------------------------------------------------

  unsigned n_field = 3;

  // Vector containers to store which boundary conditions we are applying to
  // each edge. (outlined above)
  Vector<Vector<Vector<unsigned>>> pinned_u_dofs(4, Vector<Vector<unsigned>>(3,Vector<unsigned>(free)));

  // Pin all four edges flat
  for (unsigned i = 0; i < n_field; i++)
  {
    pinned_u_dofs[0][i] = pinned_edge_yn_dof;
    pinned_u_dofs[1][i] = pinned_edge_xn_dof;
    pinned_u_dofs[2][i] = pinned_edge_yn_dof;
    pinned_u_dofs[3][i] = pinned_edge_xn_dof;
  }


  // Loop over all the boundaries in our bulk mesh
  unsigned n_bound = Bulk_mesh_pt->nboundary();
  for (unsigned b = 0; b < n_bound; b++)
  {
    // Number of elements on b
    const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);

    // Loop over the elements on boundary b
    for (unsigned e = 0; e < nb_element; e++)
    {
      // Get pointer to bulk element adjacent to b
      ELEMENT* el_pt =
	dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b, e));

      // Loop over each of the three fields
      for(unsigned i_field = 0; i_field < n_field; i_field++)
      {
	// Number of dofs we are pinning on boundary b
	const unsigned n_pinned_u_dofs = pinned_u_dofs[b][i_field].size();

	// Pin in-plane dofs (enumerated as explained above) for all nodes on
	// boundary b. Here we're applying homogeneous BCs so all pinned values
	// are simply set to zero using Parameters::get_null_fct.
	for (unsigned k = 0; k < n_pinned_u_dofs; k++)
	{
		unsigned k_type = pinned_u_dofs[b][i_field][k];
		std::cout << "On boundary " << b
				<< " pinning deflection " << i_field
				<< " type " << k_type << std::endl;
		el_pt->set_boundary_condition(
			i_field, k_type, b, Parameters::get_null_fct);
	} // end for loop over types that need to be pinned [k]
      } // end for loop over displacements [i_field]
    } // end for loop over elements on b [e]
  } // end for loop over boundaries [b]

} // end set bc



//==start_of_damped_solve ================================================
/// Used damped solves to get close to a steady solution, when close
/// enough, attempt a steady solve. If that fails, be stricter about the
/// meaning of "close enough" and repeat until a steady solve succeeds.
///
/// Expects:
///   dt_supplied_guess -- a guess for a good timestep size
///             epsilon -- an 'error tolerance' for the timestepper to limit
///                        the size of a damped step
/// Accepts:
///        doc_unsteady -- do we document every unsteady (nonphysical) solve
///   begin_with_steady -- do we begin by trying a steady solve
///
/// Returns:
///   A tuple containing suggestions to use as the inputs in the next
///   damped_solve. First a double suggestion for dt which is the timestep
///   of the first successful solve -- if the size of deformation is roughly
///   the same, it should work fine. Secondly a bool suggestion for
///   begin_with_steady which is true if this damped_solve took <2 timesteps.
//========================================================================
template<class ELEMENT>
std::tuple<double, bool> UnstructuredKSProblem<ELEMENT>::damped_solve(
  const double& dt_supplied_guess,
  const double& epsilon,
  const bool& doc_unsteady,
  const bool& begin_with_steady)
{
  // We are unsteady until a steady solve succeeds
  bool steady = false;
  // Max residual of the steady problem before we attempt a steady solve
  double sufficiently_small = 1.0e-2;
  // Timestep size
  double dt = dt_supplied_guess;
  // Try steady indicates when we should attempt a steady solve
  bool try_steady = begin_with_steady;
  // Value to be returned for initial guess for next damped_solve dt
  double dt_return_guess = dt;
  // Only set dt_initial_guess once, after the first successful solve
  bool dt_return_guess_is_unset = true;

  // If we are documenting the damped stage, create an initial state before any
  // unsteady solves have been done.
  if (doc_unsteady)
  {
    doc_solution(false);
  }
  else
  {
    Doc_unsteady_info.number()++;
  }

  while (!steady)
  {
    //------------------------------------------------------------------------
    // If we are supposed to try a steady solve, do it
    // [zdec] maybe move this if() into a helper function that is called as
    // appropriate
    if (try_steady)
    {
      oomph_info << "ATTEMPT A STEADY SOLVE" << std::endl;
      // Store the dofs before a steady solve so that they can be put back in
      // case it fails
      store_current_dof_values();
      // Get the max residual in case we need it to adjust sufficiently_small
      DoubleVector res;
      get_residuals(res);
      double max_steady_residual = res.max();
      try
      {
	// <<< Solve >>>
        steady_newton_solve();
	// If that worked, we have achieved steady state.
	// Celebrate and take note
	oomph_info << "\nHOORAY\n"
		   << "Steady solve was successful, damped solve complete\n"
		   << std::endl;
        steady = true;
	// If we are documenting the unsteady states, add the final solution to
	// the unsteady solution outputs
        if (doc_unsteady)
        {
          doc_solution(false);
        }
      }
      // If the steady solve fails, we need to tidy up before carrying on
      catch (OomphLibError& error)
      {
	// If our tolerance to attempt a steady solve is smaller than the
	// tolerance, then this implies that the initial residual was within
	// tolerance and we still got an error! SHIT!!
        if (sufficiently_small < newton_solver_tolerance())
        {
          oomph_info << "\nUH OH\n"
		     << "\"sufficiently small\" is now " << sufficiently_small
                     << " which is smaller than the solver tolerance.\n"
                     << "Our initial residual was smaller than the tolerance"
		     << " and we still got an error which is bloody stupid.\n"
		     << "Giving up on damped solves..." << std::endl;
          throw error;
        } // End of if tolerance is to small
        else
        {
          oomph_info << "\nNOT STEADY ENOUGH.\n"
		     << "\"sufficiently_small\" is insufficiently small so we\n"
		     << "are decreasing it from " << sufficiently_small
		     << " to " << max_steady_residual / 2.0 << " from now on.\n"
		     << "Returning to damped solves...\n"
		     << std::endl;
	  // Decrease the threshold for attempting steady solves as this one
	  // didn't work
          sufficiently_small = max_steady_residual / 2.0;
	  // Go back to the state we were in before attempting the steady solve
          restore_dof_values();
          for (unsigned i = 0; i < ntime_stepper(); i++)
          {
            time_stepper_pt(i)->undo_make_steady();
          }
	  // Stop trying steady solves
	  try_steady = false;
	  // Keep calm and carry on
          error.disable_error_message();

        } // End of else tolerance is not too small
      } // End of catch error
    } // End of if try_steady

    //------------------------------------------------------------------------
    // Try get us close to a steady solution by solving the damped version of
    // the equations. When it is time to try a steady solve, break this loop.
    while(!try_steady)
    {
      //----------------------------------------------------------------------
      // Begin by doing a damped solve
      oomph_info << "NEW DAMPED PSEUDO-TIME STEP WITH: dt = "
		 << dt << std::endl;
      double dt_next = adaptive_unsteady_newton_solve(dt, epsilon);
      dt = dt_next;

      // If we haven't set the initial guess for the next damped solve dt, then
      // set it now. It should be the recommended timestep after the first
      // successful solve. Assuming the following damped solve will start in a
      // roughly similar state to this one, this is appropriate.
      if (dt_return_guess_is_unset)
      {
	dt_return_guess = dt_next;
	dt_return_guess_is_unset = false;
      }

      // If we are documenting the unsteady solutions then do so, else just
      // just increase the unsteady step counter to keep count of damped steps
      if (doc_unsteady)
      {
	doc_solution(false);
      }
      else
      {
	Doc_unsteady_info.number()++;
      }

      //------------------------------------------------------------------------
      // Check how close we are to a steady solution by getting the steady
      // max residual, if it is sufficiently small, try a steady solve.
      // If that doesn't work, restrict what it means to be "sufficiently small"
      // and return to unsteady. We repeat this until the steady solve works,
      // or, we give up.

      // First set the timesteppers to steady
      for (unsigned i = 0; i < ntime_stepper(); i++)
      {
	time_stepper_pt(i)->make_steady();
      }

      // Then get the residual
      DoubleVector res;
      get_residuals(res);
      double max_steady_residual = res.max();
      oomph_info << std::endl
		 << "The max steady residual is " << max_steady_residual
		 << std::endl;

      // If it is "sufficiently small" then try a steady solve
      try_steady = max_steady_residual < sufficiently_small;

      // Reset time steppers
      for (unsigned i = 0; i < ntime_stepper(); i++)
      {
	time_stepper_pt(i)->undo_make_steady();
      }
    } // End of while(!try_steady)
  } // End of while(!steady)

  //----------------------------------------------------------------------
  // Steady state achieved, return suggestions for next damped solve
  // -- A good timestep to use would be the first successful one
  // -- Next time start with a steady solve if it took one or zero damped steps
  bool next_time_begin_with_steady = Doc_unsteady_info.number() - 1 < 2;

  return std::make_tuple(dt_return_guess, next_time_begin_with_steady);
}



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::doc_solution(bool steady,
                                                   const std::string& comment)
{
  if (MPI_Helpers::communicator_pt()->my_rank() == 0)
  {
    ofstream some_file;
    char filename[100];

    // Dump to a restart file if we are steady
    if (steady)
    {
      // Write restart file
      sprintf(filename,
              "%s/restart%i.dat",
              Parameters::output_dir.c_str(),
              Doc_steady_info.number());
      some_file.open(filename);
      dump_it(some_file);
      some_file.close();
    }

    unsigned npts;
    // // Number of plot points for coarse output
    // npts = 2;
    // if(steady)
    //  {
    //   sprintf(filename, "%s/coarse_soln_%i.dat",
    //			 Parameters::output_dir.c_str(),
    //			 Doc_steady_info.number());
    //  }
    // else
    //  {
    //   sprintf(filename, "%s/coarse_soln_%i_%i.dat",
    //			 Parameters::output_dir.c_str(),
    //			 Doc_steady_info.number(),
    //			 Doc_unsteady_info.number());
    //  }
    // some_file.open(filename);
    // Bulk_mesh_pt->output(some_file,npts);
    // some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
    //			 << comment << "\"\n";
    // some_file.close();

    // Number of plot points for fine outpout
    npts = 5;
    if (steady)
    {
      sprintf(filename,
              "%s/soln_%i.dat",
              Parameters::output_dir.c_str(),
              Doc_steady_info.number());
    }
    else
    {
      sprintf(filename,
              "%s/soln_%i_%i.dat",
              Parameters::output_dir.c_str(),
              Doc_steady_info.number(),
              Doc_unsteady_info.number());
    }
    some_file.open(filename);
    Bulk_mesh_pt->output(some_file, npts);
    some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" << comment << "\"\n";
    some_file.close();
    ParaviewHelper::write_pvd_information(
      *(Parameters::pvd_stream_pt), filename, time());


    // Write the pressure, degree of swelling and
    // deflection to the trace file
    //-----------------------------------------------------
    // Get the centre deflection first
    Vector<double> origin(2, 0.0), s(2, 0.0);
    GeomObject* centre_element_pt;
    Vector<Vector<double>> u_centre(3, Vector<double>(6, 0.0));
    MeshAsGeomObject mesh_as_geom(Bulk_mesh_pt);

    mesh_as_geom.locate_zeta(origin, centre_element_pt, s);
    if (centre_element_pt != NULL)
    {
      dynamic_cast<ELEMENT*>(centre_element_pt)
	->interpolated_koiter_steigmann_disp(s, u_centre);
    }

    Trace_file_nondim << Doc_steady_info.number() << " "
                      << Doc_unsteady_info.number() << " " << time() << " "
                      << Parameters::P_mag << " "
                      << Parameters::C_swell_data_pt->value(0) << " "
					<< u_centre[0][0] << " "
					<< u_centre[1][0] << " "
					<< u_centre[2][0] << " "
                      << endl;

    Trace_file_dim << Doc_steady_info.number() << " "
                   << Doc_unsteady_info.number() << " " << time() << " "
                   << Parameters::P_mag * Parameters::P_dim << " "
                   << Parameters::C_swell_data_pt->value(0) << " "
                   << u_centre[0][0] * Parameters::L_dim << " "
                   << u_centre[1][0] * Parameters::L_dim << " "
                   << u_centre[2][0] * Parameters::L_dim << " "
			 << endl;

    //   // Are we storing XZ profile slices?
    //   if(steady && STORE_XZ_PROFILES)
    //   {
    //     // Number of evenly spaced slices between y=0(included) and y=1(not
    //     included) unsigned n_yplanes=10;
    //     // Number of evenly spaced points per slice
    //     unsigned n_ppoint=1001;
    //     Vector<double> ppoint(2,0.0), ppoint_s(2,0.0), z(8,0.0);
    //     GeomObject* ppoint_element_pt;
    //     double h = Parameters::L1/(n_ppoint-1.0);

    //     // Loop over the y-values for the xz slice planes
    //     for(unsigned i=0; i<n_yplanes; i++)
    //     {
    //	double yi = 0.5 * (double)i / (double)10.0;

    //	// Store a profile deflection slice through y=yi
    //	//-----------------------------------------------------
    //	sprintf(filename,
    //		"%s/xz_profile_y%.2f_%i.dat",
    //		Parameters::output_dir.c_str(),
    //		yi,
    //		Doc_steady_info.number());
    //	some_file.open(filename);

    //	ppoint[0]=-Parameters::L1/2.0;
    //	ppoint[1]=yi;
    //	for(unsigned j=0; j<=n_ppoint; j++)
    //	{
    //		for(unsigned i=0; i<n_element; i++)
    //		{
    //			dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i))
    //				->locate_zeta(ppoint, ppoint_element_pt, ppoint_s);
    //			if(ppoint_element_pt!=NULL)
    //			{
    //				z = dynamic_cast<ELEMENT*>(ppoint_element_pt)
    //		->interpolated_u_koiter_steigmann(ppoint_s);
    //				some_file << ppoint[0]+z[6] << " " << ppoint[1]+z[7] << " " <<
    // z[0] << std::endl;					break;
    //			}
    //		}
    //		ppoint[0]+=h;
    //	}
    //	some_file.close();
    //     }
    //   }
  }

  // Increment the doc_info numbers
  if (steady)
  {
    Doc_steady_info.number()++;
    Doc_unsteady_info.number() = 0;
  }
  else
  {
    Doc_unsteady_info.number()++;
  }
} // end of doc

//============start_of_delete_flux_elements==============================
/// Delete Poisson Flux Elements and wipe the surface mesh
//=======================================================================
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::delete_traction_elements(
  Mesh* const& surface_mesh_pt)
{
  // How many surface elements are in the surface mesh
  unsigned n_element = surface_mesh_pt->nelement();

  // Loop over the surface elements
  for (unsigned e = 0; e < n_element; e++)
  {
    // Kill surface element
    delete surface_mesh_pt->element_pt(e);
  }

  // Wipe the mesh
  surface_mesh_pt->flush_element_and_node_storage();

} // end of delete_flux_elements


//=======start_of_main========================================
/// Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char** argv)
{
  // Initialise oomph-lib's MPI
  MPI_Helpers::init(argc, argv);

  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  // Define possible command line arguments and parse the ones that
  // were actually specified
  // Directory for solution
  Parameters::output_dir = "RESLT";
  CommandLineArgs::specify_command_line_flag("--dir", &Parameters::output_dir);

  string restart_file_name;
  // Are we restarting from a dumped file?
  CommandLineArgs::specify_command_line_flag("--restart", &restart_file_name);

  // Channel length Ratio
  CommandLineArgs::specify_command_line_flag("--L", &Parameters::L1);

  // Poisson Ratio
  CommandLineArgs::specify_command_line_flag("--nu", &Parameters::Nu);

  // Damping coefficient
  CommandLineArgs::specify_command_line_flag("--mu", &Parameters::Mu);

  // Prestrain on the membrane
  CommandLineArgs::specify_command_line_flag("--prestrain",
					     &Parameters::Isotropic_prestrain);

  // Prestrain max to end simultation
  double prestrain_max = 0.03;
  CommandLineArgs::specify_command_line_flag("--prestrain_max", &prestrain_max);

  // Prestrain increment
  double prestrain_inc = 0.001;
  CommandLineArgs::specify_command_line_flag("--prestrain_inc", &prestrain_inc);

  // Increment size for Pressure
  double p_inc_dim = 500.0;
  CommandLineArgs::specify_command_line_flag("--p_inc", &p_inc_dim);

  // Incrememnt size for swelling
  double p_max_dim = 2000.0;
  CommandLineArgs::specify_command_line_flag("--p_max", &p_max_dim);

  // Maximum degree of swelling
  CommandLineArgs::specify_command_line_flag("--c_max", &Parameters::C_swell_max);

  // Minimum degree of deswelling
  CommandLineArgs::specify_command_line_flag("--c_min", &Parameters::C_swell_min);

  // Incrememnt size for swelling
  CommandLineArgs::specify_command_line_flag("--c_inc", &Parameters::C_swell_inc);

  // Element Area (no larger element than)
  CommandLineArgs::specify_command_line_flag("--element_area",
                                             &Parameters::element_area);

  // Use linear stress rather than MR
  CommandLineArgs::specify_command_line_flag("--use_linear_stress");

  // Use the finite difference jacobian
  CommandLineArgs::specify_command_line_flag("--use_fd_jacobian");

  // Pin u_alpha everywhere
  CommandLineArgs::specify_command_line_flag("--pininplane");

  // Run the debug block at the end of main
  CommandLineArgs::specify_command_line_flag("--debug");


  // Parse command line
  CommandLineArgs::parse_and_assign();

  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();

  // How many nodes do we want to manually place along the boundaries
  // (roughly length*width/thickness)
  Parameters::n_long_edge_nodes = ceil(0.0 * Parameters::L1 / Parameters::Thickness) + 2;
  Parameters::n_short_edge_nodes = ceil(0.0 * Parameters::L2 / Parameters::Thickness) + 2;

  // What are the dimensionless versions of the pressure parameters
  Parameters::P_inc = p_inc_dim / Parameters::P_dim;
  Parameters::P_max = p_max_dim / Parameters::P_dim;


  //===========================================================================
  // CREATE THE PROBLEM
  UnstructuredKSProblem<KoiterSteigmannC1CurvableBellElement> problem;

  // Set up some problem paramters
  problem.newton_solver_tolerance() = 1e-9;
  problem.max_residuals() = 1e4;
  problem.max_newton_iterations() = 10;
  problem.target_error_safety_factor() = 0.5;
  //problem.enable_xz_profiles();

  // Restart the problem if a restart file was passed
  if (CommandLineArgs::command_line_flag_has_been_set("--restart"))
  {
    oomph_info << "We are restarting from: " << restart_file_name << std::endl;
    ifstream restart_stream;
    restart_stream.open(restart_file_name);
    problem.restart(restart_stream);
    restart_stream.close();
    problem.Doc_steady_info.number() -= 1;
  }

  // Set up timestepping parameters and history values
  double dt = problem.next_dt();
  bool try_steady = false;
  double epsilon = 1.0e-3;
  problem.assign_initial_values_impulsive();
  problem.initialise_dt(dt);

  // Open the pvd file that tracks all the solutions
  ofstream pvd_stream;
  pvd_stream.open(Parameters::output_dir + "/steady_solns.pvd");
  ParaviewHelper::write_pvd_header(pvd_stream);
  Parameters::pvd_stream_pt = &pvd_stream;

  // Perform an initial solve to get the state before any computations
  oomph_info
    << "=================================================================\n"
    << "=================================================================\n"
    << "DO AN INITIAL STATE SOLVE\n"
    << "=================================================================\n"
    << "=================================================================\n"
    << std::endl;
  problem.steady_newton_solve(); // SOLVE


  // Pre-inflate the membrane
  oomph_info
    << "=================================================================\n"
    << "=================================================================\n"
    << "INFLATION STAGE\n"
    << "=================================================================\n"
    << "=================================================================\n"
    << std::endl;
  // If we aren't restarting, do the inflation step
  if (!CommandLineArgs::command_line_flag_has_been_set("--restart"))
  {
    while(Parameters::P_mag<Parameters::P_max)
    {
      // INFLATION
      Parameters::P_mag = Parameters::P_max;
      // <<< Solve >>>
      std::tie(dt, try_steady) =
      problem.damped_solve(dt, epsilon, false, try_steady);
      // Document steady
      problem.doc_solution(true); // AND DOCUMENT
    }
  }


  // Pre-inflate the membrane
  oomph_info
    << "=================================================================\n"
    << "=================================================================\n"
    << "PRESTRAINING STAGE\n"
    << "=================================================================\n"
    << "=================================================================\n"
    << std::endl;
  // First release prestrain
  Parameters::Isotropic_prestrain = 0.0;
  // <<< Solve >>>
  std::tie(dt, try_steady) =
    problem.damped_solve(dt, epsilon, false, try_steady);
  // Document steady
  problem.doc_solution(true);
  // Then incrememnt prestrain to our limit
  while(Parameters::Isotropic_prestrain < prestrain_max)
  {
    // Prestrain
    Parameters::Isotropic_prestrain += prestrain_inc;
    // <<< Solve >>>
    std::tie(dt, try_steady) =
      problem.damped_solve(dt, epsilon, false, try_steady);
    // Document steady
    problem.doc_solution(true); // AND DOCUMENT
  }



  // Close the pvd file
  ParaviewHelper::write_pvd_footer(pvd_stream);
  pvd_stream.close();

  // Print success
  oomph_info << "Exiting Normally\n";

  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize();
} // End of main
