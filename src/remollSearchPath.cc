#include "remollSearchPath.hh"

#include <G4GenericMessenger.hh>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/param.h>
#include <limits.h>
#include <unistd.h>

remollSearchPath* remollSearchPath::gInstance = 0;
remollSearchPath* remollSearchPath::GetInstance() {
  if (gInstance == 0) {
    gInstance = new remollSearchPath();
  }
  return gInstance;
}

remollSearchPath::remollSearchPath()
: fPrefix("")
{
  // Determine prefix
  Initialize();

  // Create messenger
  fSearchPathMessenger = new G4GenericMessenger(this,"/remoll/searchpath/","Remoll search path properties");
}

remollSearchPath::~remollSearchPath()
{
  delete fSearchPathMessenger;
}

void remollSearchPath::Initialize()
{
  // Find executable prefix
  char buffer[MAXPATHLEN];
  #ifdef __APPLE__
  uint32_t len = sizeof(buffer);
  if (_NSGetExecutablePath(buffer, &len) == 0) {
    fPrefix = std::string(buffer);
  } else G4cerr << "ERROR - Could not read executable path, assuming \"./\"" << G4endl;
  #else // Linux
  ssize_t len = readlink("/proc/self/exe", buffer, MAXPATHLEN-1);
  if (len != -1) {
    buffer[len] = '\0';
    fPrefix = std::string(buffer);
  } else G4cerr << "ERROR - Could not read executable path, assuming \"./\"" << G4endl;
  #endif
  fPrefix = fPrefix.substr(0,fPrefix.rfind('/')); // strip `remoll`
  fPrefix = fPrefix.substr(0,fPrefix.rfind('/')); // strip `bin` or `build`
  G4cout << "Using remoll resources under " << fPrefix << G4endl;

  // Set possible subdirs
  fSubDirs.push_back("");
  fSubDirs.push_back("/share");
  fSubDirs.push_back("/share/remoll");
}

std::string remollSearchPath::Find(std::string file)
{
  G4cout << "Searching " << file << G4endl;
  for (size_t i = 0; i < fSubDirs.size(); i++) {
    G4cout << "in " << fSubDirs.at(i) << G4endl;
    struct stat info;
    std::string path = fPrefix + fSubDirs.at(i);
    if (stat(path.c_str(),&info) == 0) {
      if (S_ISREG(info.st_mode)) {
        G4cout << "Using " << path << G4endl;
        return path;
      }
    }
  }
  return file;
}
