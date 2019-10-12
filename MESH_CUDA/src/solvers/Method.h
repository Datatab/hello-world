/*
 * Method.h
 *
 *  Created on: Oct 11, 2019
 *      Author: v1
 */

#ifndef METHOD_H_
#define METHOD_H_

class Method {
public:
	virtual void init(char* xmlFileName) = 0;
	virtual void run() = 0;
	virtual void done() = 0;
};

#endif /* METHOD_H_ */
