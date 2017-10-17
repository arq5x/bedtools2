/*
 * ContextSubtract.h
 *
 *  Created on: Feb 19, 2015
 *      Author: nek3d
 */

#ifndef CONTEXTSUBTRACT_H_
#define CONTEXTSUBTRACT_H_


#include "ContextIntersect.h"

class ContextSubtract : public ContextIntersect {
public:
	ContextSubtract();
	virtual ~ContextSubtract();
	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);
    virtual bool hasIntersectMethods() const { return true; }
    virtual bool isValidState();
    void setSubtractFraction(float fraction) { _fractionalSubtract = fraction; }
    float getSubtractFraction() const { return _fractionalSubtract; }
    bool getRemoveAll() const { return _removeAll; }
    bool getRemoveSum() const { return _removeSum; }


private:
    bool handle_A();
    bool handle_N();

    float _fractionalSubtract;
    bool _removeAll;
    bool _removeSum;

};

#endif /* CONTEXTSUBTRACT_H_ */
