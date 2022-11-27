/*
 * bastructs.h
 *
 * Copyright (C) 2013 Louis-Francois Handfield
 * e-mail: lfhandfield@gmail.com
 *
 * This program is free software; upon notification by email to the licensor
 * of the licencee identity and nature of use, the licencee can redistribute
 * this program and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2
 * of the License, or (at the licencee option) any later version. As such,
 * no further notifications are required.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 *
 * This file contains template specialization and non-template
 *
 */

// optionnal external dependency, which uses a modified file "GSLfunc.hpp"

#ifndef _defined_Cstruct
#define _defined_Cstruct

#define GNU_SCIENTIFIC_LIBRARY
#ifndef LFH_HAS_RUNNINGTIME_STATISTICS
#define LFH_HAS_RUNNINGTIME_STATISTICS true
#endif

#include "core.h"
#include "bastructs.h"

using namespace std;
namespace LFHPrimitive{

	char* cloneString(const char* const what);
// Anything can be int double char*, and optionally have an integer index
// The type itself holds the size of an unit within the 0xFF byte

// structure store a list of strings.
// structure maintain a scope for substring searches


} // end of namespace LFHPrimitive

#endif


