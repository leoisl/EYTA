//
// Created by Leandro Ishi Soares de Lima on 07/02/17.
//

#include "EnhanceTranscriptomeTool.h"

void EnhanceTranscriptomeTool::checkParameters() {
  //check the transcriptome files
  string transcriptsFile = getInput()->getStr(STR_PAR_TRANSCRIPTOME);
  checkIfFileExists(transcriptsFile);
}