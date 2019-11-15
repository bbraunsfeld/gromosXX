/**
 * @file qm_atom.h
 * QM atom structure for QM/MM
 */

#ifndef INCLUDED_QM_ATOM_H
#define	INCLUDED_QM_ATOM_H

namespace interaction {
  /**
   * @struct QM_Atom
   * A structure holding information on QM atoms
   */
  struct QM_Atom {
    /**
     * Constructor
     * @param index - the index of the atom in topology
     * @param atomic_number - the atomic number of the atom
     * @param pos - 3D vector of atomic position
     * @param linked - the index of atom it is linked to - default is 0 (none)
     */
    QM_Atom(unsigned index
         , unsigned atomic_number = 1
         , math::Vec pos = {0.0,0.0,0.0}
         ) : index(index)
           , atomic_number(atomic_number)
           , pos(pos)
           , force(0.0)
           , qm_charge(0.0)
    {}
    /**
     * Copy constructor
     */
    QM_Atom(const QM_Atom & a) : index(a.index)
                               , atomic_number(a.atomic_number)
                               , pos(a.pos)
                               , force(a.force)
                               , qm_charge(a.qm_charge)
    {}

    /**
     * the index of the atom in topology - starts with 0
     */
    const unsigned index;
    
    /**
     * the atomic number of the atom
     */
    const unsigned atomic_number;

    /**
     * The coordinate
     */
    mutable math::Vec pos;

    /**
     * The force
     */
    mutable math::Vec force;

    /**
     * The charge calculated from QM
     */
    mutable double qm_charge;

    /**
     * less-than comparison operator
     */
    bool operator<(const QM_Atom & a) const {
      return index < a.index;
    }
  };
}

#endif	/* QM_ATOM_H */

