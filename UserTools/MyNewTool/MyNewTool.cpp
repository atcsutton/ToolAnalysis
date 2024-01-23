#include "MyNewTool.h"

MyNewTool::MyNewTool():Tool(){}


bool MyNewTool::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  return true;
}


bool MyNewTool::Execute(){

  std::cout << "Hello world" << std::endl;

  return true;
}


bool MyNewTool::Finalise(){

  return true;
}
