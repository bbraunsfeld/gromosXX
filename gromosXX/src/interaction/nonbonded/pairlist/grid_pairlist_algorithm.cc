/**
 * @file grid_pairlist_algorithm.cc
 * grid pairlist algorithm
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction_types.h>
#include <math/periodicity.h>
#include <math/volume.h>

#include <interaction/nonbonded/pairlist/pairlist.h>
#include <interaction/nonbonded/interaction/storage.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

#include <interaction/nonbonded/interaction/nonbonded_term.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_term.h>

#include <interaction/nonbonded/interaction/nonbonded_innerloop.h>

#include <interaction/nonbonded/pairlist/pairlist_algorithm.h>
#include <interaction/nonbonded/pairlist/grid_pairlist_algorithm.h>

// #include <interaction/nonbonded/interaction_spec.h>
// #include <util/template_split.h>

#include <util/debug.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

interaction::Grid_Pairlist_Algorithm::Grid_Pairlist_Algorithm()
  : interaction::Pairlist_Algorithm()
{
}

int interaction::Grid_Pairlist_Algorithm::init
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 bool quiet)
{

  set_cutoff(sim.param().pairlist.cutoff_short, 
	     sim.param().pairlist.cutoff_long);

  grid_properties(topo, conf, sim);

  if (!quiet){
    std::cout << "GridPairlistAlgorithm\n"
	      << "\tcells             " 
	      << std::setw(10) << m_grid.Na
	      << std::setw(10) << m_grid.Nb
	      << std::setw(10) << m_grid.Nc << "\n"
	      << "\textenced cells    "
	      << std::setw(10) << m_grid.Na_ex
	      << std::setw(10) << m_grid.Nb_ex
	      << std::setw(10) << m_grid.Nc_ex << "\n"
	      << "\tcell size         "
	      << std::setw(10) << m_grid.a
	      << std::setw(10) << m_grid.b
	      << std::setw(10) << m_grid.c << "\n";

    const int Ncell = m_grid.Na * m_grid.Nb * m_grid.Nc;
    const int N = m_grid.Na * m_grid.Nb;

    const double P = topo.num_chargegroups();
    const double V = math::volume(conf.current().box, conf.boundary_type);
    // const double Vcell = m_grid.a * m_grid.b * m_grid.c;

    const double Vcut = 4.0 / 3.0 * math::Pi * m_cutoff_long_2 * m_cutoff_long;
    const double Pcut = P * Vcut / V;
    
    const double Pcell = P / Ncell;
    const double Player = P / m_grid.Nc;
    
    std::cout << "\tparticles / cell    " << std::setw(10) << Pcell << "\n"
	      << "\tparticles / layer   " << std::setw(10) << Player << "\n";

    // mask size:
    int Nmask = 0;
    for(int z=0; z<=m_grid.mask_z; ++z){
      for(unsigned int y=0; y<m_grid.mask[z].size(); y+=2){
	Nmask += m_grid.mask[z][y+1] - m_grid.mask[z][y] - 1;
      }
    }
    
    std::cout << "\tcells in mask       " << std::setw(10) << Nmask << "\n"
	      << "\tpairs               " << std::setw(10) << 0.5 * P * P << "\n"
	      << "\tpairs (grid)        " << std::setw(10) << P * Nmask * Pcell << "\n"
	      << "\tpairs (cutoff)      " << std::setw(10) << 0.5 * P * Pcut << "\n";

    // just for fun, already try this
    prepare_grid(topo, conf, sim);
    // collapse_grid();

    int occupied = 0;
    for(int z=0; z<m_grid.Nc; ++z){
      for(int i=0; i<N; ++i){
	if (m_grid.count[z][i])
	  ++occupied;
      }
    }

    std::cout << "\toccupied            " << std::setw(10) << occupied << "\n"
	      << "\tparticle / occ cell " << std::setw(10) << double(P) / occupied << "\n";
    

    // print_mask();

    std::cout << "END\n";
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// put chargegroups into box and on the grid
// the grid will be sparse...
////////////////////////////////////////////////////////////////////////////////

/**
 * put the chargegroups into the box and on the grid
 */
void interaction::Grid_Pairlist_Algorithm::prepare_grid
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  const int Nab = m_grid.Na * m_grid.Nb;
  
  // this is enough for small cells (when only about one particle fits in anyway)
  const int space = m_grid.Pcell * 3;

  m_grid.p_cell.resize(m_grid.Nc);
  m_grid.cell_start.resize(m_grid.Nc);
  m_grid.count.resize(m_grid.Nc);
  m_grid.cell_index.resize(m_grid.Nc);
  
  for(int z=0; z<m_grid.Nc; ++z){
    m_grid.p_cell[z].resize(Nab * space);
    m_grid.count[z].assign(Nab, 0);
    m_grid.cell_start[z].resize(Nab+1);
    
    int b = 0;
    for(int i=0; i<Nab; ++i, b+=space)
      m_grid.cell_start[z][i] = b;
  }
  
  math::Periodicity<math::rectangular> periodicity(conf.current().box);

  math::VArray &pos = conf.current().pos;
  math::Vec v, v_box, trans;

  topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin(),
    cg_to = topo.chargegroup_end();

  // solute chargegroups...
  unsigned int i = 0;
  for( ; i < topo.num_solute_chargegroups(); ++cg_it, ++i){
    // cog
    cg_it.cog(pos, v);

    // gather on first atom...
    v_box = v;
    periodicity.put_into_box(v_box);
    trans = v_box - v;

    // now grid the cg
    const int x = int((v_box(0) + 0.5 * conf.current().box(0)(0)) / m_grid.a);
    const int y = int((v_box(1) + 0.5 * conf.current().box(1)(1)) / m_grid.b);
    const int z = int((v_box(2) + 0.5 * conf.current().box(2)(2)) / m_grid.c);
  
    const int c = y * m_grid.Na + x;

    m_grid.p_cell[z][m_grid.cell_start[z][c] + m_grid.count[z][c]] = Grid::Particle(i, v_box);
    ++m_grid.count[z][c];
    
    // atoms in a chargegroup
    topology::Atom_Iterator at_it = cg_it.begin(),
      at_to = cg_it.end();
    for( ; at_it != at_to; ++at_it){
      pos(*at_it) += trans;

    } // loop over atoms
  } // loop over solute cg's

  // solvent chargegroups
  for( ; cg_it != cg_to; ++cg_it, ++i){

    // cog is first atom
    v = pos(**cg_it);
    v_box = v;
    periodicity.put_into_box(v_box);
    trans = v_box - v;
    
    // now grid the cg
    const int x = int((v_box(0) + 0.5 * conf.current().box(0)(0)) / m_grid.a);
    const int y = int((v_box(1) + 0.5 * conf.current().box(1)(1)) / m_grid.b);
    const int z = int((v_box(2) + 0.5 * conf.current().box(2)(2)) / m_grid.c);
    
    const int c = y * m_grid.Na + x;
    
    m_grid.p_cell[z][m_grid.cell_start[z][c] + m_grid.count[z][c]] = Grid::Particle(i, v_box);
    ++m_grid.count[z][c];

    // loop over the atoms
    topology::Atom_Iterator at_it = cg_it.begin(),
      at_to = cg_it.end();
    for( ; at_it != at_to; ++at_it){
      pos(*at_it) += trans;
    } // atoms
  } // solvent cg's

  // check that there was enough space
  for(int z=0; z < m_grid.Nc; ++z){
    for(int i=0; i<Nab; ++i)
      if (m_grid.count[z][i] > space){
	io::messages.add("Not enough space to put chargegroups into cells!",
			 "Grid Pairlist Algorithm",
			 io::message::critical);
      }
  }
}

void interaction::Grid_Pairlist_Algorithm::collapse_grid()
{
  const int N = m_grid.Na * m_grid.Nb;

  for(int z=0; z<m_grid.Nc; ++z){
    int p_cell_index = m_grid.count[z][0];
    int cell_index_ex = (m_grid.Nb_ex - m_grid.Nb) / 2 * m_grid.Na_ex +
      (m_grid.Na_ex - m_grid.Na) / 2;

    m_grid.cell_index[z].clear();
    
    // need the index for all particles...
    for(int j=0; j<m_grid.count[z][0]; ++j)
      m_grid.cell_index[z].push_back(cell_index_ex);
    
    for(int i=1; i < N; ++i){
      const int start = m_grid.cell_start[z][i];
      m_grid.cell_start[z][i] = p_cell_index;

      ++cell_index_ex;
      if (i % m_grid.Na == 0)
	cell_index_ex += m_grid.Na_ex - m_grid.Na;

      for(int j=0; j < m_grid.count[z][i]; ++j, ++p_cell_index){

	m_grid.p_cell[z][p_cell_index] = m_grid.p_cell[z][start + j];
	m_grid.cell_index[z].push_back(cell_index_ex);

      } // loop over particles in cell

    } // loop over cells in plane

    // sentinel...
    m_grid.cell_start[z][N] = p_cell_index;

  } // loop over planes (z)
  
}

/**
 * check the grid and print it...
 */
void interaction::Grid_Pairlist_Algorithm::print_grid()
{
  std::cout << "THE GRID\n========\n";

  const int N = m_grid.Na * m_grid.Nb;
  int num_P = 0;
  int errors = 0;

  const double ha = 0.5 * m_grid.Na * m_grid.a;
  const double hb = 0.5 * m_grid.Nb * m_grid.b;
  const double hc = 0.5 * m_grid.Nc * m_grid.c;

  for(int z=0; z<m_grid.Nc; ++z){
    
    // std::cout << "plane " << z << ":\n";
    int ci = 0;
    
    for(int i=0; i<N; ++i){
      if (m_grid.count[z][i] != m_grid.cell_start[z][i+1] - m_grid.cell_start[z][i]){
	std::cout << "grid error: count and cell_start don't match!" << std::endl;
	++errors;
      }
      
      if (m_grid.count[z][i] != 0){
	
	
	for(int n=0; n<m_grid.count[z][i]; ++n, ++ci){

	  // look up extended index
	  const int ind = m_grid.cell_index[z][ci];
	  const int indy = ind / m_grid.Na_ex - (m_grid.Nb_ex - m_grid.Nb) / 2;
	  const int indx = ind % m_grid.Na_ex - (m_grid.Na_ex - m_grid.Na) / 2;
	  if (indy * m_grid.Na + indx != i){
	    std::cout << "grid error: extended index and index don't match" << std::endl;
	    ++errors;
	  }

	  Grid::Particle const & p = m_grid.p_cell[z][m_grid.cell_start[z][i] + n];
	  
	  if (p.shift_index != 0){
	    std::cout << "grid error: shift index not zero in central box" << std::endl;
	    ++errors;
	  }
	  
	  if (p.z + hc < z * m_grid.c || p.z + hc > (z+1) * m_grid.c ||
	      p.x + ha < indx * m_grid.a || p.x + ha > (indx+1) * m_grid.a ||
	      p.y + hb < indy * m_grid.b || p.y + hb > (indy+1) * m_grid.b){
	    std::cout << "grid error: particle cog not inside cell" << std::endl;
	    // std::cout << "x=" << p.x << " y=" << p.y << " z=" << p.z << "\n"
	    // << "indx*a=" << indx*m_grid.a << " indy*b=" << indy*m_grid.b
	    // << " z*c=" << z*m_grid.c << "\n";
	    ++errors;
	  }
	  
	  std::cout << "cell " << i << " (" << ind << " = [" << indx << ", " << indy << "]) ["
		    << m_grid.count[z][i] << "] ";

	  std::cout << p.i << "\n";
	  
	  ++num_P;
	}
	
	std::cout << "\n";
      }
    }
  }
  
  std::cout << "particles on grid: " << num_P << "\n";
  
  if (errors == 0){
    std::cout << "no errors detected\n";
  }
  else
    std::cout << errors << " errors detected!\n";
}

/**
 * calculate grid properties,
 * put chargegroups on grid,
 * including the center of geometries
 */
void interaction::Grid_Pairlist_Algorithm::prepare
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  DEBUG(7, "standard pairlist algorithm : prepare");
  
  set_cutoff(sim.param().pairlist.cutoff_short, 
	     sim.param().pairlist.cutoff_long);
  
  // first put the chargegroups into the box
  // _prepare_cog(conf, topo);

  grid_properties(topo, conf, sim);

  prepare_grid(topo, conf, sim);

  collapse_grid();

  // print_grid();
  
}

void interaction::Grid_Pairlist_Algorithm::grid_properties
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  const double s = sim.param().pairlist.grid_size;
  
  m_grid.Na = int(rint(conf.current().box(0)(0) / s));
  m_grid.Nb = int(rint(conf.current().box(1)(1) / s));
  m_grid.Nc = int(rint(conf.current().box(2)(2) / s));
  
  m_grid.a = conf.current().box(0)(0) / m_grid.Na;
  m_grid.b = conf.current().box(1)(1) / m_grid.Nb;
  m_grid.c = conf.current().box(2)(2) / m_grid.Nc;  

  const int Ncell = m_grid.Na * m_grid.Nb * m_grid.Nc;
  
  const double P = topo.num_chargegroups();
  // const double V = math::volume(conf.current().box, conf.boundary_type);
  // const double Vcell = m_grid.a * m_grid.b * m_grid.c;

  m_grid.Pcell = int(P / Ncell) + 1;

  m_grid.Na_ex = int(m_cutoff_long / m_grid.a);
  if (m_cutoff_long / m_grid.a > math::epsilon) ++m_grid.Na_ex;
  m_grid.Na_ex *= 2;
  m_grid.Na_ex += m_grid.Na;

  m_grid.Nb_ex = int(m_cutoff_long / m_grid.b);
  if (m_cutoff_long / m_grid.b > math::epsilon) ++m_grid.Nb_ex;
  m_grid.Nb_ex *= 2;
  m_grid.Nb_ex += m_grid.Nb;

  m_grid.Nc_ex = int(m_cutoff_long / m_grid.c);
  if (m_cutoff_long / m_grid.c > math::epsilon) ++m_grid.Nc_ex;
  m_grid.Nc_ex *= 2;
  m_grid.Nc_ex += m_grid.Nc;

  calculate_mask();

  // and the shift vectors
  m_shift_vector.clear();

  for(int z=0; z<2; ++z){
    for(int y=-1; y<2; ++y){
      for(int x=-1; x<2; ++x){
	m_shift_vector.push_back(x * conf.current().box(0) + y * conf.current().box(1) + z * conf.current().box(2));
      }
    }
  }

}

/**
 * calculate the mask
 */
void interaction::Grid_Pairlist_Algorithm::calculate_mask()
{
  const double c2 = m_grid.c * m_grid.c;
  const double b2 = m_grid.b * m_grid.b;
  // const double a2 = m_grid.a * m_grid.a;

  m_grid.mask_z = int(m_cutoff_long / m_grid.c);
  if (m_cutoff_long / m_grid.c > math::epsilon) ++m_grid.mask_z;
  
  m_grid.mask.resize(m_grid.mask_z+1);
  

  double z_dist, y_dist;

  // special case of 0 plane
  {
    m_grid.mask[0].clear();
    int mask_y = int(sqrt(m_cutoff_long_2) / m_grid.b);
    if (sqrt(m_cutoff_long_2) / m_grid.b > math::epsilon) ++mask_y;
    
    for(int y=0; y <= mask_y; ++y){
      const int row = y * m_grid.Na_ex;

      if (y>1) y_dist = (y - 1) * (y - 1) * b2;
      else y_dist = 0.0;

      int mask_x = int(sqrt(m_cutoff_long_2 - y_dist) / m_grid.a);
      if (sqrt(m_cutoff_long_2 - y_dist) / m_grid.a > math::epsilon) ++mask_x;
      // begin 'till one past end
      m_grid.mask[0].push_back(row - mask_x);
      m_grid.mask[0].push_back(row + mask_x + 1);
    }
    // don't do self interaction over mask...
    m_grid.mask[0][0] = 1;
  }

  for(int z=1; z<=m_grid.mask_z; ++z){
    m_grid.mask[z].clear();

    if (z>1) z_dist = (z - 1) * (z - 1) * c2;
    else z_dist = 0.0;
    
    int mask_y = int(sqrt(m_cutoff_long_2 - z_dist) / m_grid.b);
    if (sqrt(m_cutoff_long_2 - z_dist) / m_grid.b > math::epsilon) ++mask_y;
    
    for(int y=mask_y; y>=0; --y){
      const int row = -y * m_grid.Na_ex;

      if (y>1) y_dist = (y - 1) * (y - 1) * b2;
      else y_dist = 0.0;

      int mask_x = int(sqrt(m_cutoff_long_2 - y_dist - z_dist) / m_grid.a);
      if (sqrt(m_cutoff_long_2 - y_dist - z_dist) / m_grid.a > math::epsilon) ++mask_x;

      // begin 'till one past end
      assert(m_grid.mask.size() > unsigned(z));
      m_grid.mask[z].push_back(row - mask_x);
      m_grid.mask[z].push_back(row + mask_x + 1);
    }
  }

  for(int z=1; z <= m_grid.mask_z; ++z){
    int row = 0;
    for(int y = m_grid.mask[z].size() - 4; y >= 0; y-=2){
      row += 2 * m_grid.Na_ex;
      m_grid.mask[z].push_back(row + m_grid.mask[z][y]);
      m_grid.mask[z].push_back(row + m_grid.mask[z][y+1]);
    }
  }

}

void interaction::Grid_Pairlist_Algorithm::print_mask()
{
  std::cout << "\tmask\n";
  for(int z=0; z <= m_grid.mask_z; ++z){
    std::cout << "\n\tplane " << z << ":\n\t";

    for(unsigned int y=0; y < m_grid.mask[z].size(); y+=2){
      
      assert(m_grid.mask.size() > unsigned(z));
      assert(m_grid.mask[z].size() > y+1);

      std::cout << std::setw(5) << m_grid.mask[z][y] << " -> " 
		<< std::setw(5) << m_grid.mask[z][y+1] << "\t";

      if ((y + 2) % 10 == 0) std::cout << "\n\t";

    }
  }
  std::cout << "\n";
}

////////////////////////////////////////////////////////////////////////////////
// prepare a plane
////////////////////////////////////////////////////////////////////////////////
void interaction::Grid_Pairlist_Algorithm::prepare_plane
(
 int z,
 std::vector<Grid::Particle> & p_plane, 
 std::vector<int> & cell_start
 )
{
  int z_shift = 0;
  if (z >= m_grid.Nc){
    z_shift = 9;
    z -= m_grid.Nc;
  }

  // reserve enough space for anything...
  p_plane.resize(m_grid.p_cell[z].size() * 4);
  cell_start.resize(m_grid.Na_ex * m_grid.Nb_ex + 1);

  const int a_ex = (m_grid.Na_ex - m_grid.Na) / 2;
  const int b_ex = (m_grid.Nb_ex - m_grid.Nb) / 2;
  
  // index into the (newly constructed) plane
  int j = 0;
  int cj = 0;
  int cs = 0;
  
  /*
    --------------------    1 = i
    |j |            |  |    2 = i_ex_to
    |  |            |  |    3 = i_ex
    --------------------    4 = i_to
    |  |            |  |
    |  |            |  |    first copy i_ex -> i_to to j
    |  |1 2      3 4|  |    then  copy i -> i_to
    |  |            |  |    then  copy i -> i_ex_to
    --------------------
    |  |            |  |    then add Na to i, i_ex_to, i_ex, i_to
    |  |            |  |
    --------------------
  */


  // the upper extended area
  int i = (m_grid.Nb - b_ex) * m_grid.Na;
  int i_to = i + m_grid.Na;
  int i_ex = i_to - a_ex;
  int i_ex_to = i + a_ex;

  for(int e=0; e<b_ex; ++e){

    const int pi = m_grid.cell_start[z][i];
    const int pi_to = m_grid.cell_start[z][i_to];
    const int pi_ex = m_grid.cell_start[z][i_ex];
    const int pi_ex_to = m_grid.cell_start[z][i_ex_to];

    // std::cout << "i=" << i << " i_to=" << i_to << " i_ex=" << i_ex 
    // << " i_ex_to=" << i_ex_to << "\n";
    // std::cout << "pi=" << pi << " pi_to=" << pi_to << " pi_ex=" << pi_ex 
    // << " pi_ex_to=" << pi_ex_to << "\n";
    
    // shift the particles
    // std::cerr << "particles" << std::endl;
    // std::cout << "from " << i_ex << " to " << i_to << ":\n";
    for(int p=pi_ex; p < pi_to; ++p, ++j){
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift, m_shift_vector);
      // std::cout << "shift " << m_grid.p_cell[z][p].i << " to cell " << j << "\n";
    }
    // std::cout << "from " << i << " to " << i_to << ":\n";
    for(int p=pi; p < pi_to; ++p, ++j){
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift + 1, m_shift_vector);
      // std::cout << "shift " << m_grid.p_cell[z][p].i << " to cell " << j << "\n";
    }
    // std::cout << "from " << i << " to " << i_ex_to << ":\n";
    for(int p=pi; p < pi_ex_to; ++p, ++j){
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift + 2, m_shift_vector);
      // std::cout << "shift " << m_grid.p_cell[z][p].i << " to cell " << j << "\n";
    }

    // adapt the cell_start
    // std::cerr << "cell start" << std::endl;
    for(int c=i_ex; c < i_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i_ex] + cs;
    cs += m_grid.cell_start[z][i_to] - m_grid.cell_start[z][i_ex];
    for(int c=i; c < i_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i] + cs;
    cs += m_grid.cell_start[z][i_to] - m_grid.cell_start[z][i];
    for(int c=i; c < i_ex_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i] + cs;
    cs += m_grid.cell_start[z][i_ex_to] - m_grid.cell_start[z][i];
    // std::cerr << "cj after upper (" << e << ") = " << cj << std::endl;

    i += m_grid.Na;
    i_to += m_grid.Na;
    i_ex += m_grid.Na;
    i_ex_to += m_grid.Na;
  }

  // the center area
  // std::cerr << "center" << std::endl;
  
  i = 0;
  i_to = i + m_grid.Na;
  i_ex = i_to - a_ex;
  i_ex_to = i + a_ex;

  for(int e=0; e<m_grid.Nb; ++e){

    const int pi = m_grid.cell_start[z][i];
    const int pi_to = m_grid.cell_start[z][i_to];
    const int pi_ex = m_grid.cell_start[z][i_ex];
    const int pi_ex_to = m_grid.cell_start[z][i_ex_to];
    
    // shift particles
    for(int p=pi_ex; p < pi_to; ++p, ++j)
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift + 3, m_shift_vector);
    for(int p=pi; p < pi_to; ++p, ++j)
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift + 4, m_shift_vector);
    for(int p=pi; p < pi_ex_to; ++p, ++j)
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift + 5, m_shift_vector);

    // adapt the cell_start
    for(int c=i_ex; c < i_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i_ex] + cs;
    cs += m_grid.cell_start[z][i_to] - m_grid.cell_start[z][i_ex];
    for(int c=i; c < i_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i] + cs;
    cs += m_grid.cell_start[z][i_to] - m_grid.cell_start[z][i];
    for(int c=i; c < i_ex_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i] + cs;
    cs += m_grid.cell_start[z][i_ex_to] - m_grid.cell_start[z][i];
    // std::cerr << "cj after upper (" << e << ") = " << cj << std::endl;

    i += m_grid.Na;
    i_to += m_grid.Na;
    i_ex += m_grid.Na;
    i_ex_to += m_grid.Na;
  }

  // and the final part
  i = 0;
  i_to = i + m_grid.Na;
  i_ex = i_to - a_ex;
  i_ex_to = i + a_ex;

  for(int e=0; e<b_ex; ++e){

    const int pi = m_grid.cell_start[z][i];
    const int pi_to = m_grid.cell_start[z][i_to];
    const int pi_ex = m_grid.cell_start[z][i_ex];
    const int pi_ex_to = m_grid.cell_start[z][i_ex_to];
    
    // shift particles
    for(int p=pi_ex; p < pi_to; ++p, ++j)
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift + 6, m_shift_vector);
    for(int p=pi; p < pi_to; ++p, ++j)
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift + 7, m_shift_vector);
    for(int p=pi; p < pi_ex_to; ++p, ++j)
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift + 8, m_shift_vector);

    // adapt the cell_start
    for(int c=i_ex; c < i_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i_ex] + cs;
    cs += m_grid.cell_start[z][i_to] - m_grid.cell_start[z][i_ex];
    for(int c=i; c < i_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i] + cs;
    cs += m_grid.cell_start[z][i_to] - m_grid.cell_start[z][i];
    for(int c=i; c < i_ex_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i] + cs;
    cs += m_grid.cell_start[z][i_ex_to] - m_grid.cell_start[z][i];
    // std::cerr << "cj after upper (" << e << ") = " << cj << std::endl;

    i += m_grid.Na;
    i_to += m_grid.Na;
    i_ex += m_grid.Na;
    i_ex_to += m_grid.Na;
  }

  // sentinel
  // std::cerr << "cj " << cj << std::endl;
  cell_start[cj] = cs;

}

/**
 * check the grid and print it...
 */
void interaction::Grid_Pairlist_Algorithm::print_plane
(
 int z,
 std::vector<Grid::Particle> & p_plane, 
 std::vector<int> & cell_start
)
{
  std::cout << "PLANE " << std::setw(3) << z << "\n=========\n";

  std::cout.precision(3);
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  
  const int N = m_grid.Na_ex * m_grid.Nb_ex;
  int num_P = 0;

  // const double ha = 0.5 * m_grid.Na * m_grid.a;
  // const double hb = 0.5 * m_grid.Nb * m_grid.b;
  // const double hc = 0.5 * m_grid.Nc * m_grid.c;

  int ci = 0;
    
  for(int i=0; i<N; ++i){
    
    if (cell_start[i+1] - cell_start[i] != 0){
	
      std::cout << "cell " << i << " [" << cell_start[i+1] - cell_start[i] << "] = ";
      
      for(int n = cell_start[i]; n < cell_start[i+1]; ++n){
	
	Grid::Particle const & p = p_plane[n];
	
	std::cout << p.i << " (" << p.x << " | " << p.y << " | " << p.z << " : " << p.shift_index << ") ";
	++num_P;
      }
      
      std::cout << "\n";
	++ci;
    }
  }
  
  std::cout << "particles on plane: " << num_P << "\n";
}


////////////////////////////////////////////////////////////////////////////////
// pairlist update
////////////////////////////////////////////////////////////////////////////////

void interaction::Grid_Pairlist_Algorithm::update
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 interaction::Storage & storage,
 interaction::Pairlist & pairlist,
 unsigned int begin,
 unsigned int end,
 unsigned int stride
 )
{
  Nonbonded_Innerloop innerloop(*m_param);
  innerloop.init(sim);

  // empty the pairlist
  for(unsigned int i=0; i<topo.num_atoms(); ++i)
    pairlist[i].clear();
  
  const double update_start = util::now();

  std::vector<Grid::Particle> p_plane;
  std::vector<int> cell_start;

  const int c_ex = (m_grid.Nc_ex - m_grid.Nc) / 2;
  const int N = m_grid.Na * m_grid.Nb;
  const int num_solute_cg = topo.num_solute_chargegroups();

  for(int z=begin; z < m_grid.Nc + c_ex; z+=stride){

    prepare_plane(z, p_plane, cell_start);
    // print_plane(z, p_plane, cell_start);

    // std::cerr << "plane " << z << " printed..." << std::endl;
    
    if (z < m_grid.Nc){
      // do self interaction
      const int i_first = m_grid.cell_start[z][0];
      const int i_last  = m_grid.cell_start[z][N];

      int base = -1;
      int start = -1;

      for(int i=i_first; i < i_last; ++i){

	if (m_grid.p_cell[z][i].i < num_solute_cg){ // self interaction
	  DEBUG(8, "self " << m_grid.p_cell[z][i].i);
	  for(int a1 = topo.chargegroups()[m_grid.p_cell[z][i].i],
		a_to = topo.chargegroups()[m_grid.p_cell[z][i].i + 1];
	      a1 < a_to; ++a1){
	    for(int a2 = a1+1; a2 < a_to; ++a2){
	      if (excluded_solute_pair(topo, a1, a2))
		continue;
	      pairlist[a1].push_back(a2);
	    }
	  }
	}
	
	if (base != m_grid.cell_index[z][i]){
	  base = m_grid.cell_index[z][i];
	  start = i;
	}
	else{ // more than 1 cg in the cell
	  if (m_grid.p_cell[z][i].i < num_solute_cg){ // solute - ?
	    for(int j=start; j<i; ++j){
	      DEBUG(8, "intra cell " << m_grid.p_cell[z][i].i << " - " 
		    << m_grid.p_cell[z][j].i);
	      
	      if (m_grid.p_cell[z][j].i < num_solute_cg){ // solute - solute

		const int ii = (m_grid.p_cell[z][i].i < m_grid.p_cell[z][j].i) ? 
		  m_grid.p_cell[z][i].i : m_grid.p_cell[z][j].i;
		const int jj = (m_grid.p_cell[z][i].i < m_grid.p_cell[z][j].i) ? 
		  m_grid.p_cell[z][j].i : m_grid.p_cell[z][i].i;

		for(int a1 = topo.chargegroups()[ii],
		      a_to = topo.chargegroups()[ii + 1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[jj],
			a2_to = topo.chargegroups()[jj + 1];
		      a2 < a2_to; ++a2){
		    
		    if (excluded_solute_pair(topo, a1, a2))
		      continue;
		    pairlist[a1].push_back(a2);
		  }
		}
	      }
	      else{ // solute - solvent
		for(int a1 = topo.chargegroups()[m_grid.p_cell[z][i].i],
		      a_to = topo.chargegroups()[m_grid.p_cell[z][i].i+1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[m_grid.p_cell[z][j].i],
			a2_to = topo.chargegroups()[m_grid.p_cell[z][j].i + 1];
		      a2 < a2_to; ++a2){
		    pairlist[a1].push_back(a2);
		  }
		}

	      }
	    }
	  }
	  else{ // solvent - ?
	    for(int j=start; j<i; ++j){
	      DEBUG(8, "intra cell " << m_grid.p_cell[z][i].i << " - " 
		    << m_grid.p_cell[z][j].i);
	      
	      for(int a1 = topo.chargegroups()[m_grid.p_cell[z][i].i],
		    a_to = topo.chargegroups()[m_grid.p_cell[z][i].i+1];
		  a1 < a_to; ++a1){
		for(int a2 = topo.chargegroups()[m_grid.p_cell[z][j].i],
		      a2_to = topo.chargegroups()[m_grid.p_cell[z][j].i + 1];
		    a2 < a2_to; ++a2){
		  pairlist[a2].push_back(a1);
		}
	      }
	    }
	  }
	}
      } // loop over plane
    } // self interaction

    ////////////////////////////////////////////////////////////////////////////////
    // 
    // INTER CELL INTERACTIONS
    //
    ////////////////////////////////////////////////////////////////////////////////

    // loop over all mask levels (inside box)
    for(int mask_z=0; mask_z < int(m_grid.mask.size()); ++mask_z){
      if (z - mask_z < 0) break;
      if (z - mask_z >= m_grid.Nc) continue;
      
      // std::cerr << "z=" << z << " mask_z=" << mask_z
      // << " i_level=" << z - mask_z << std::endl;
      
      const int i_level = z - mask_z;
      
      // std::cout << "mask level " << mask_z << std::endl;
      // std::cerr << "mask level " << mask_z << std::endl;

      const int i_first = m_grid.cell_start[i_level][0];
      const int i_last  = m_grid.cell_start[i_level][N];

      // std::cout << "particles " << i_last - i_first << std::endl;
      // std::cerr << "particles " << i_last - i_first << std::endl;

      // loop over all particles in this level
      for(int i=i_first; i < i_last; ++i){
	
	assert(m_grid.cell_index[i_level].size() > unsigned(i));
	const int base = m_grid.cell_index[i_level][i];
	
	if (m_grid.p_cell[i_level][i].i < num_solute_cg){ // solute - ?
	  
	  assert(m_grid.mask.size() > unsigned(mask_z));
	  for(unsigned int m=0; m<m_grid.mask[mask_z].size(); m+=2){
	    
	    assert(cell_start.size() > unsigned(base + m_grid.mask[mask_z][m]));
	  
	    const int first = cell_start[base + m_grid.mask[mask_z][m]];
	    const int last  = cell_start[base + m_grid.mask[mask_z][m+1]];
	    
	    for(int j=first; j<last; ++j){
	      
	      const double d2 = 
		(m_grid.p_cell[i_level][i].x - p_plane[j].x) *
		(m_grid.p_cell[i_level][i].x - p_plane[j].x) +
		(m_grid.p_cell[i_level][i].y - p_plane[j].y) *
		(m_grid.p_cell[i_level][i].y - p_plane[j].y) +
		(m_grid.p_cell[i_level][i].z - p_plane[j].z) *
		(m_grid.p_cell[i_level][i].z - p_plane[j].z);

	      if (d2 > m_cutoff_long_2) continue;

	      // no self interaction here...
	      DEBUG(8, "inter cell: " << m_grid.p_cell[i_level][i].i << " - " << p_plane[j].i);
	      assert(m_grid.p_cell[i_level][i].i != p_plane[j].i);

	      if (p_plane[j].i < num_solute_cg){ 
		// solute - solute

		if (d2 > m_cutoff_short_2){

		  for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
			a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		      a1 < a_to; ++a1){
		    for(int a2 = topo.chargegroups()[p_plane[j].i],
			  a2_to = topo.chargegroups()[p_plane[j].i + 1];
			a2 < a2_to; ++a2){
		      
		      // std::cout << "a1 = " << math::v2s(conf.current().pos(a1)) << "\n"
		      // << "a2 = " << math::v2s(conf.current().pos(a2)) << "\n"
		      // << "shift = " << math::v2s(m_shift_vector[p_plane[j].shift_index])
		      // << "\n";
		      
		      innerloop.lj_crf_innerloop_shift(topo, conf, a1, a2, storage,
						       m_shift_vector[p_plane[j].shift_index]);
		    }
		  }
		}
		else{
		  
		  const int ii = (m_grid.p_cell[i_level][i].i < p_plane[j].i) ? 
		    m_grid.p_cell[i_level][i].i : p_plane[j].i;
		  const int jj = (m_grid.p_cell[i_level][i].i < p_plane[j].i) ? 
		    p_plane[j].i : m_grid.p_cell[i_level][i].i;
		  
		  DEBUG(8, "rewritten: " << ii
			<< " - " << jj);
		  
		  for(int a1 = topo.chargegroups()[ii],
			a_to = topo.chargegroups()[ii + 1];
		      a1 < a_to; ++a1){
		    for(int a2 = topo.chargegroups()[jj],
			  a2_to = topo.chargegroups()[jj + 1];
			a2 < a2_to; ++a2){
		      
		      // std::cout << "ii=" << ii << " jj=" << jj
		      // << " a1=" << a1 << " a2=" << a2 << std::endl;
		      
		      if (excluded_solute_pair(topo, a1, a2))
			continue;
		      pairlist[a1].push_back(a2);
		    }
		  }
		}
	      }
	      else{ // solute - solvent
		if (d2 > m_cutoff_short_2){

		  for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
			a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		      a1 < a_to; ++a1){
		    for(int a2 = topo.chargegroups()[p_plane[j].i],
			  a2_to = topo.chargegroups()[p_plane[j].i + 1];
			a2 < a2_to; ++a2){
		      
		      innerloop.lj_crf_innerloop_shift(topo, conf, a1, a2, storage,
						       m_shift_vector[p_plane[j].shift_index]);
		    }
		  }
		}
		else{

		  for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
			a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		      a1 < a_to; ++a1){
		    for(int a2 = topo.chargegroups()[p_plane[j].i],
			  a2_to = topo.chargegroups()[p_plane[j].i + 1];
			a2 < a2_to; ++a2){
		      pairlist[a1].push_back(a2);
		    }
		  }
		}
		
	      }
	      
	    } // j in mask row
	  }  // mask
	} // solute - ?
	else{ // solvent - ?
	  assert(m_grid.mask.size() > unsigned(mask_z));
	  for(unsigned int m=0; m<m_grid.mask[mask_z].size(); m+=2){
	    
	    assert(cell_start.size() > unsigned(base + m_grid.mask[mask_z][m]));
	  
	    const int first = cell_start[base + m_grid.mask[mask_z][m]];
	    const int last  = cell_start[base + m_grid.mask[mask_z][m+1]];
	    
	    for(int j=first; j<last; ++j){

	      const double d2 = 
		(m_grid.p_cell[i_level][i].x - p_plane[j].x) *
		(m_grid.p_cell[i_level][i].x - p_plane[j].x) +
		(m_grid.p_cell[i_level][i].y - p_plane[j].y) *
		(m_grid.p_cell[i_level][i].y - p_plane[j].y) +
		(m_grid.p_cell[i_level][i].z - p_plane[j].z) *
		(m_grid.p_cell[i_level][i].z - p_plane[j].z);

	      if (d2 > m_cutoff_long_2) continue;

	      DEBUG(8, "inter (solvent - ?): " << m_grid.p_cell[i_level][i].i << " - " << p_plane[j].i);
	      DEBUG(8, "i_level=" << i_level << " d2=" << d2);
	      
	      if (d2 > m_cutoff_short_2){
		for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
		      a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[p_plane[j].i],
			a2_to = topo.chargegroups()[p_plane[j].i + 1];
		      a2 < a2_to; ++a2){

		    // maybe for fast solvent loop, check cg of j...
		    innerloop.lj_crf_innerloop_shift(topo, conf, a1, a2, storage,
						     m_shift_vector[p_plane[j].shift_index]);
		  }
		}
	      }
	      else{
		for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
		      a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[p_plane[j].i],
			a2_to = topo.chargegroups()[p_plane[j].i + 1];
		      a2 < a2_to; ++a2){
		    pairlist[a2].push_back(a1);
		  }
		}
	      }
	      
	    }
	    
	  }

	} // solvent - ?

	// std::cout << "\n";
	
      } // particle in the current mask plane

    } // mask levels

  } // the extended planes

  m_timing += util::now() - update_start;
}

bool interaction::Grid_Pairlist_Algorithm
::excluded_solute_pair(topology::Topology & topo,
		       unsigned int i, unsigned int j)
{
  assert(i<j);
  
  std::set<int>::const_reverse_iterator
    e = topo.all_exclusion(i).rbegin(),
    e_to = topo.all_exclusion(i).rend();

  for( ; e != e_to; ++e){
    if (j > unsigned(*e)) break;
    if (j == unsigned(*e)){
      DEBUG(11, "\texcluded");
      return true;
    }
      
  }
  DEBUG(12, "\tnot excluded");
  return false;
}