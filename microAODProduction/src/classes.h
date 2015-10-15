/*!
 * \file classes.h
 * \brief Definition of PatProduction classes for which ROOT dictionary should be made.
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Maria Teresa Grippo (University of Siena, INFN Pisa)
 * \date 2014-04-27 created
 *
 * Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
 *                Maria Teresa Grippo <grippomariateresa@gmail.com>
 *
 * This file is part of X->HH->bbTauTau.
 *
 * X->HH->bbTauTau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * X->HH->bbTauTau is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "HHbbTauTau/PatProduction/interface/PatVertex.h"

namespace {
struct dictionary {
    pat::Vertex vertex;
    std::vector<pat::Vertex>::const_iterator v_p_v_ci;
    edm::Wrapper<std::vector<pat::Vertex> >  w_v_p_v;
    std::vector<float> v_f;
    std::vector<math::XYZVector> v_3d;
};
}
