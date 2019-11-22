/**
 * @file qm_link.cc
 * Implements methods of QM_Link
 */
#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>

#include <math/boundary_implementation.h>
#include <math/periodicity.h>

#include <util/debug.h>


#include <interaction/qmmm/qm_atom.h>
#include <interaction/qmmm/mm_atom.h>
#include <interaction/qmmm/qm_link.h>
#include <interaction/qmmm/qm_zone.h>

//#include <util/template_split.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

void interaction::QM_Link::update_cap_position(const math::Vec& qm_pos, const math::Vec& mm_pos, const double cap_length) const
  {
  math::Vec r_qmmm = mm_pos - qm_pos;
  DEBUG(15, "QM pos: " << math::v2s(qm_pos));
  DEBUG(15, "MM pos: " << math::v2s(mm_pos));
  DEBUG(15, "Link vector: " << math::v2s(r_qmmm));
  DEBUG(15, "Cap length: " << cap_length);
  
  this->pos = qm_pos + cap_length * (r_qmmm / math::abs(r_qmmm));
  DEBUG(15, "Capping atom position: " << math::v2s(this->pos));
}

void interaction::QM_Link::distribute_force(interaction::QM_Zone& qm_zone) const
  {
  DEBUG(15, "Redistributing force of capping atom");
  DEBUG(15, "QM-MM link: " << this->qm_index << " - " << this->mm_index);
  std::set<QM_Atom>::iterator qm_it = qm_zone.qm.find(this->qm_index);
  std::set<MM_Atom>::iterator mm_it = qm_zone.mm.find(this->mm_index);
  assert(qm_it != qm_zone.qm.end());
  assert(mm_it != qm_zone.mm.end());
  math::Vec& qm_pos = qm_it->pos,
             mm_pos = mm_it->pos,
             qm_force = qm_it->force,
             mm_force = mm_it->force;
  math::Vec r_mm_qm = mm_pos - qm_pos;

  DEBUG(15, "QM pos: " << math::v2s(qm_pos));
  DEBUG(15, "QM force: " << math::v2s(qm_force));
  DEBUG(15, "MM pos: " << math::v2s(mm_pos));
  DEBUG(15, "MM force: " << math::v2s(mm_force));
  DEBUG(15, "r_mm_qm: " << math::v2s(r_mm_qm));
  

  double d_l_qm = math::abs(this->pos - qm_pos),
         d2_mm_qm = math::abs2(r_mm_qm),
         d_mm_qm = sqrt(d2_mm_qm),
         link_fraction = d_l_qm / d_mm_qm;
  
  math::Vec lr_r_d2 = link_fraction * r_mm_qm / d2_mm_qm;
  DEBUG(15, "Force on capping atom: " << math::v2s(this->force));
  for (unsigned i = 0; i < 3; ++i) {
    DEBUG(15, "Coordinate " << i);
    math::Vec a_qm(0.0);
    math::Vec a_mm(0.0);
    a_qm[i] = 1 - link_fraction;
    DEBUG(15, "a_qm: " << math::v2s(a_qm));
    a_mm[i] = link_fraction;
    DEBUG(15, "a_mm: " << math::v2s(a_mm));
    math::Vec b = (mm_pos(i) - qm_pos(i)) * lr_r_d2;
    DEBUG(15, "a_qm + b: " << math::v2s(a_qm + b));
    DEBUG(15, "a_mm - b: " << math::v2s(a_mm - b));
    qm_force[i] += math::dot(this->force, a_qm + b);
    mm_force[i] += math::dot(this->force, a_mm - b);
    DEBUG(15, "Force[" << i << "] distributed to QM atom: " << math::dot(this->force, a_qm + b));
    DEBUG(15, "Force[" << i << "] distributed to MM atom: " << math::dot(this->force, a_mm - b));
  }
  // Zero capping atom force, as it has been distributed
  this->force = 0.0;
}