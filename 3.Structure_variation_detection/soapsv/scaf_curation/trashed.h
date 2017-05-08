/*
 * trashed.h
 *
 *  Created on: Feb 24, 2010
 *      Author: Aqua
 */

#ifndef TRASHED_H_
#define TRASHED_H_

/*
Trashes

if(marker_pending && (line_vec[0].find("@SQ") != string::npos))
{
	marker_pending = false;
	GOSSIP("SQ specific");
	if(!__GetSVInfomation(line_vec, param_ptr->sv))
	{
		DBG(string("@SQ header missing in file: ") + param_ptr->filename);
		continue;
	}
	sv_left_inc = (param_ptr->sv.start_pos + 1);
	sv_right_exc = (sv_left_inc + param_ptr->sv.genotype.size());
}

*/

#endif /* TRASHED_H_ */
