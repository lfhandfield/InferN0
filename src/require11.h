/*

 Copyright (C) 2013 Louis-Francois Handfield
 All rights reserved.

 */
#ifndef _defined_LFHRequire11
#define _defined_LFHRequire11

#include "bastructs.h"
#include "primitive.h"

#define CRT_THREAD(ObJeCt,FuNcTiOnNaMe,...) std::thread(&decltype(ObJeCt)::FuNcTiOnNaMe, &ObJeCt,##__VA_ARGS__)

namespace LFHPrimitive{
template<class SCOPE, class FUNCTION, class OUTPUT, typename... ARGUMENTS>
class FunctionCall{
public:

};
} // end of namespace

#include "require11.hpp"

#endif
