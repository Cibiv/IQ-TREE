/*
 * mhrw.h
 *
 *  Created on: Sept 1, 2018
 *      Author: Lukasz
 */
#pragma once

#include <string>
#include <set>
#include <vector>
#include "srw.h"

class MHRW : public SRW {
public:
    
protected:
    void transition() override;
};

