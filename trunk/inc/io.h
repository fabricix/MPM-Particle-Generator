/*
 * io.h
 *
 *  Created on: Oct 6, 2018
 *      Author: fabricio
 */

#ifndef INC_IO_H_
#define INC_IO_H_

namespace IO
{
    void ReadInputFile(int argc, char **argv);
    void WriteOutputFile();
    std::string pathGet();
}
#endif /* INC_IO_H_ */
