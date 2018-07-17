#ifndef __REMOLLSEARCHPATH_HH
#define __REMOLLSEARCHPATH_HH 1

#include <vector>
#include <string>

class G4GenericMessenger;

class remollSearchPath {
  private:
    // Singleton pointer
    static remollSearchPath* gInstance;
    // Private constructor
    remollSearchPath();

  public:
    // Public destructor
    virtual ~remollSearchPath();
    // Static instance getter
    static remollSearchPath* GetInstance();

  private:
    // Initialization
    void Initialize();

  private:
    // Prefix derived from executable
    std::string fPrefix;
    // Ordered list of additional paths to check under fPrefix
    std::vector<std::string> fSubDirs;

    // Messenger
    G4GenericMessenger* fSearchPathMessenger;

  public:
    std::string Find(std::string file);
    const char* Find(const char* file) {
      std::string path = GetInstance()->Find(std::string(file));
      return path.c_str();
    };
};

#endif // __REMOLLSEARCHPATH_HH
