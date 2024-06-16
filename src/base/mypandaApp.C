#include "mypandaApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
mypandaApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_bahavior") = false;
  return params;
}

mypandaApp::mypandaApp(InputParameters parameters) : MooseApp(parameters)
{
  mypandaApp::registerAll(_factory, _action_factory, _syntax);
}

mypandaApp::~mypandaApp() {}

void
mypandaApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<mypandaApp>(f, af, s);
  Registry::registerObjectsTo(f, {"mypandaApp"});
  Registry::registerActionsTo(af, {"mypandaApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
mypandaApp::registerApps()
{
  registerApp(mypandaApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
mypandaApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  mypandaApp::registerAll(f, af, s);
}
extern "C" void
mypandaApp__registerApps()
{
  mypandaApp::registerApps();
}
