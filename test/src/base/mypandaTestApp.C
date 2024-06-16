//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "mypandaTestApp.h"
#include "mypandaApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
mypandaTestApp::validParams()
{
  InputParameters params = mypandaApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

mypandaTestApp::mypandaTestApp(InputParameters parameters) : MooseApp(parameters)
{
  mypandaTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

mypandaTestApp::~mypandaTestApp() {}

void
mypandaTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  mypandaApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"mypandaTestApp"});
    Registry::registerActionsTo(af, {"mypandaTestApp"});
  }
}

void
mypandaTestApp::registerApps()
{
  registerApp(mypandaApp);
  registerApp(mypandaTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
mypandaTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  mypandaTestApp::registerAll(f, af, s);
}
extern "C" void
mypandaTestApp__registerApps()
{
  mypandaTestApp::registerApps();
}
