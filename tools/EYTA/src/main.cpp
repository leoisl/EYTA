/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

// We include the header file for the tool
#include "build_dbg.hpp"
#include "map_reads.hpp"
#include "EnhanceTranscriptome.h"
#include "EnhanceTranscriptomeTool.h"
#include "global.h"
#include "Tests.h"
#include "Utils.h"
/********************************************************************************/

int main (int argc, char* argv[])
{
    try
    {
        //Build DBG
        cout << "Step 1/ Building Merged DBG..." << endl;
        build_dbg().run (argc, argv); //this call will set up graph and nodeIdToUnitigId
        cout << "Done!" << endl;
        cout.flush();

        //mapping and ULG
        cout << "Step 2/ Mapping transcripts and building ULG..." << endl;
        map_reads().run (argc, argv);
        cout << "Done!" << endl;
        cout.flush();

        //Find alternative splicing events
        cout << "Step 3/ Finding alternative splicing events..." << endl;
        EnhanceTranscriptomeTool().run (argc, argv);
        cout << "Done!" << endl;
        cout.flush();
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

