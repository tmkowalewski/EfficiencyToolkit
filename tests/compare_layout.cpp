#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TKey.h>

#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

////////////////////////////////////////////////////////////////////////////////
/// Compare all saved objects and numerical graph values in a detector directory.
void CompareDetector(TDirectory &previous, TDirectory &current)
{
   std::set<std::string> expected;
   for (auto object : *previous.GetListOfKeys()) {
      auto *key = static_cast<TKey *>(object);
      const std::string name = key->GetName();
      expected.insert(name);
      auto *oldObject = previous.Get(name.c_str());
      auto *newObject = current.Get(name.c_str());
      if (!oldObject || !newObject || oldObject->IsA() != newObject->IsA())
         throw std::runtime_error("Missing or different object: " + name);
      if (auto *oldGraph = dynamic_cast<TGraphErrors *>(oldObject)) {
         auto *newGraph = static_cast<TGraphErrors *>(newObject);
         if (oldGraph->GetN() != newGraph->GetN())
            throw std::runtime_error("Graph point count differs: " + name);
         for (int i = 0; i < oldGraph->GetN(); ++i)
            if (oldGraph->GetPointX(i) != newGraph->GetPointX(i) || oldGraph->GetPointY(i) != newGraph->GetPointY(i) ||
                oldGraph->GetErrorX(i) != newGraph->GetErrorX(i) || oldGraph->GetErrorY(i) != newGraph->GetErrorY(i))
               throw std::runtime_error("Graph values differ: " + name);
      }
      if (name == "normalization" && std::string(oldObject->GetTitle()) != newObject->GetTitle())
         throw std::runtime_error("Normalization differs");
   }
   std::set<std::string> actual;
   for (auto object : *current.GetListOfKeys())
      actual.insert(object->GetName());
   if (expected != actual)
      throw std::runtime_error("Detector object names differ");
}

////////////////////////////////////////////////////////////////////////////////
/// Check that the old grouped output and new flat output contain equivalent detectors.
int main(int argc, char **argv)
{
   if (argc != 3) {
      std::cerr << "Usage: compare_layout grouped.root flat.root\n";
      return 1;
   }
   try {
      TFile oldFile(argv[1], "READ"), newFile(argv[2], "READ");
      if (oldFile.IsZombie() || newFile.IsZombie())
         throw std::runtime_error("Cannot open comparison files");
      std::set<std::string> expected;
      for (const char *group : {"clover_cross", "clover_back", "cebr_all"}) {
         auto *directory = oldFile.GetDirectory(group);
         if (!directory || newFile.GetDirectory(group))
            throw std::runtime_error(std::string("Invalid group layout: ") + group);
         for (auto object : *directory->GetListOfKeys()) {
            const std::string name = object->GetName();
            if (!expected.insert(name).second)
               throw std::runtime_error("Duplicate detector name: " + name);
            auto *previous = directory->GetDirectory(name.c_str());
            auto *current = newFile.GetDirectory(name.c_str());
            if (!previous || !current)
               throw std::runtime_error("Missing detector: " + name);
            CompareDetector(*previous, *current);
         }
      }
      std::set<std::string> actual;
      for (auto object : *newFile.GetListOfKeys())
         actual.insert(object->GetName());
      if (expected != actual || expected.empty())
         throw std::runtime_error("Top-level detector names differ");
      std::cout << "Verified " << expected.size() << " detector directories and all graph data\n";
      return 0;
   } catch (const std::exception &error) {
      std::cerr << error.what() << '\n';
      return 1;
   }
}
