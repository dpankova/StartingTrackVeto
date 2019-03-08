#
# $Id: CMakeLists.txt 141949 2016-04-27 16:53 kjero $
#

i3_project(StartingTrackVetoLE
  PYTHON_DIR python)

i3_add_library(StartingTrackVetoLE
  private/StartingTrackVetoLE/StartingTrackVetoLE.cxx    
  private/StartingTrackVetoLE/StartingTrackVetoLEUtils.cxx    
  USE_PROJECTS icetray dataclasses recclasses phys-services
               photonics-service serialization
  USE_TOOLS boost python
)
