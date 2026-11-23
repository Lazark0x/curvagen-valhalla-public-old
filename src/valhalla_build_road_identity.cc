// Writes a tileset's road-identity sidecar (thor/road_twin_index.h) to a file, so the
// engine maps it at start (thor.roundtrip_road_identity_file) instead of scanning every
// tile. Run it on the finished pack: the file is bound to the tiles it was built from.
#include "argparse_utils.h"
#include "baldr/graphreader.h"
#include "midgard/logging.h"
#include "thor/road_twin_index.h"

#include <boost/property_tree/ptree.hpp>
#include <cxxopts.hpp>

#include <filesystem>
#include <iostream>

int main(int argc, char* argv[]) {
  const auto program = std::filesystem::path(__FILE__).stem().string();
  boost::property_tree::ptree config;
  std::string out;

  try {
    // clang-format off
    cxxopts::Options options(
      program,
      program + " " + VALHALLA_PRINT_VERSION + "\n\n"
      "writes the round-trip road-identity sidecar of the configured tileset to a file.\n"
      "The thor.roundtrip_* radii and switchback setting are read from the same config\n"
      "the engine will run with.\n\n");

    options.add_options()
      ("h,help", "Print this help message.")
      ("v,version", "Print the version of this software.")
      ("c,config", "Path to the json configuration file.", cxxopts::value<std::string>())
      ("i,inline-config", "Inline json config.", cxxopts::value<std::string>())
      ("o,output", "The file to write.", cxxopts::value<std::string>(out));
    // clang-format on

    auto result = options.parse(argc, argv);
    if (!parse_common_args(program, options, result, &config))
      return EXIT_SUCCESS;
    if (out.empty())
      throw cxxopts::exceptions::exception("--output is required\n\n" + options.help());
  } catch (cxxopts::exceptions::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (std::exception& e) {
    std::cerr << "Unable to parse command line options because: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  valhalla::baldr::GraphReader reader(config.get_child("mjolnir"));
  const double twin_r = config.get<double>("thor.roundtrip_twin_radius_m", 30.0);
  const double par_r = config.get<double>("thor.roundtrip_parallel_radius_m", 80.0);
  const bool parallels = config.get<bool>("thor.roundtrip_parallel_tier", true);
  const bool switchback = config.get<bool>("thor.roundtrip_switchback_test", true);

  const auto idx = valhalla::thor::RoadTwinIndex::Scan(reader, twin_r, par_r, parallels, switchback);
  if (!idx->Save(out, reader)) {
    LOG_ERROR("could not write " + out);
    return EXIT_FAILURE;
  }
  // Read it back through the engine's own loader, so a file that would be refused fails
  // here, at pack build, and not on a phone.
  std::string why;
  const auto back =
      valhalla::thor::RoadTwinIndex::Load(out, reader, twin_r, par_r, parallels, switchback, &why);
  if (!back || !back->same_contents(*idx)) {
    LOG_ERROR(out + " does not read back as written" + (why.empty() ? "" : ": " + why));
    return EXIT_FAILURE;
  }
  LOG_INFO("wrote " + out + " (" + std::to_string(std::filesystem::file_size(out)) + " bytes)");
  return EXIT_SUCCESS;
}
