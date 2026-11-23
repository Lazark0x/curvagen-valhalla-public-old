// The road-identity sidecar persisted at pack build (valhalla_build_road_identity) must
// answer exactly as the one scanned from the tiles, and a file the engine cannot trust
// must be refused, never used.

#include "baldr/graphreader.h"
#include "gurka.h"
#include "test.h"
#include "thor/road_twin_index.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace valhalla;
using valhalla::thor::RoadTwinIndex;

namespace {

std::string file_beside(const gurka::map& m, const std::string& name) {
  return (std::filesystem::path(m.config.get<std::string>("mjolnir.tile_dir")).parent_path() / name)
      .string();
}

std::vector<char> read_all(const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

void write_all(const std::string& path, const std::vector<char>& bytes) {
  std::ofstream(path, std::ios::binary).write(bytes.data(), bytes.size());
}

} // namespace

// ---------------------------------------------------------------------------------
// 4 m/char.
//
//   P================================================Q   "PQ", one-way
//   S================================================R   "SR", one-way, 8 m: TWIN
//
//   T------------------------------------------------U   "TU", two-way, 44 m: PARALLEL
// ---------------------------------------------------------------------------------
class RoadIdentityPersist : public ::testing::Test {
protected:
  static gurka::map map;
  static gurka::map other;
  static void SetUpTestSuite() {
    const std::string ascii_map = R"(
P                                                                                                   Q

S                                                                                                   R










T                                                                                                   U
    )";
    const gurka::ways ways = {
        {"PQ", {{"highway", "primary"}, {"oneway", "yes"}}},
        {"SR", {{"highway", "primary"}, {"oneway", "yes"}}},
        {"TU", {{"highway", "secondary"}}},
    };
    map = gurka::buildtiles(gurka::detail::map_to_coordinates(ascii_map, 4), ways, {}, {},
                            "test/data/road_identity_persist");

    const std::string other_map = R"(
A                                                   B

D                                                   C
    )";
    other = gurka::buildtiles(gurka::detail::map_to_coordinates(other_map, 4),
                              {{"AB", {{"highway", "primary"}, {"oneway", "yes"}}},
                               {"CD", {{"highway", "primary"}, {"oneway", "yes"}}}},
                              {}, {}, "test/data/road_identity_persist_other");
  }
};
gurka::map RoadIdentityPersist::map = {};
gurka::map RoadIdentityPersist::other = {};

TEST_F(RoadIdentityPersist, PersistedIndexAnswersLikeTheScannedOne) {
  auto reader = test::make_clean_graphreader(map.config.get_child("mjolnir"));
  const auto scanned = RoadTwinIndex::Scan(*reader, 30.0, 80.0, true, true);
  ASSERT_GT(scanned->stats().twin_keys, 0u) << "map has no twins: the comparison proves nothing";
  ASSERT_GT(scanned->stats().parallel_keys, 0u)
      << "map has no parallels: the comparison proves nothing";

  const std::string path = file_beside(map, "road_identity.bin");
  ASSERT_TRUE(scanned->Save(path, *reader));
  std::string why;
  const auto loaded = RoadTwinIndex::Load(path, *reader, 30.0, 80.0, true, true, &why);
  ASSERT_TRUE(loaded) << why;
  EXPECT_TRUE(loaded->mapped());
  EXPECT_TRUE(loaded->same_contents(*scanned));

  // Every query the engine makes, on every directed edge of the map.
  uint32_t with_twins = 0, with_parallels = 0;
  for (const auto& tid : reader->GetTileSet()) {
    auto tile = reader->GetGraphTile(tid);
    for (uint32_t i = 0; i < tile->header()->directededgecount(); ++i) {
      const baldr::GraphId e(tid.tileid(), tid.level(), i);
      const uint64_t canon =
          RoadTwinIndex::canonical_id(tile->directededge(e), e, reader->GetOpposingEdgeId(e));
      std::vector<uint64_t> st, lt, sp, lp;
      scanned->append_twins(canon, st);
      loaded->append_twins(canon, lt);
      scanned->append_parallels(canon, sp);
      loaded->append_parallels(canon, lp);
      EXPECT_EQ(st, lt) << "twins of " << canon;
      EXPECT_EQ(sp, lp) << "parallels of " << canon;
      EXPECT_EQ(scanned->has_twins(canon), loaded->has_twins(canon)) << canon;
      EXPECT_EQ(scanned->maybe_has_twins(canon), loaded->maybe_has_twins(canon)) << canon;
      with_twins += !st.empty();
      with_parallels += !sp.empty();
    }
  }
  EXPECT_GT(with_twins, 0u);
  EXPECT_GT(with_parallels, 0u);
}

TEST_F(RoadIdentityPersist, RefusesAFileItCannotTrust) {
  auto reader = test::make_clean_graphreader(map.config.get_child("mjolnir"));
  const std::string path = file_beside(map, "road_identity_refuse.bin");
  ASSERT_TRUE(RoadTwinIndex::Scan(*reader, 30.0, 80.0, true, true)->Save(path, *reader));
  ASSERT_TRUE(RoadTwinIndex::Load(path, *reader, 30.0, 80.0, true, true));

  std::string why;
  EXPECT_FALSE(RoadTwinIndex::Load(path, *reader, 25.0, 80.0, true, true, &why)) << "twin radius";
  EXPECT_FALSE(RoadTwinIndex::Load(path, *reader, 30.0, 70.0, true, true, &why)) << "parallel radius";
  EXPECT_FALSE(RoadTwinIndex::Load(path, *reader, 30.0, 80.0, false, true, &why)) << "parallels";
  EXPECT_FALSE(RoadTwinIndex::Load(path, *reader, 30.0, 80.0, true, false, &why)) << "switchback";
  EXPECT_FALSE(RoadTwinIndex::Load(path + ".missing", *reader, 30.0, 80.0, true, true, &why));

  auto other_reader = test::make_clean_graphreader(other.config.get_child("mjolnir"));
  EXPECT_FALSE(RoadTwinIndex::Load(path, *other_reader, 30.0, 80.0, true, true, &why))
      << "a file from other tiles was accepted";
  EXPECT_NE(why.find("other tiles"), std::string::npos) << why;

  const auto bytes = read_all(path);
  const std::string bad = file_beside(map, "road_identity_bad.bin");
  write_all(bad, std::vector<char>(bytes.begin(), bytes.end() - 4));
  EXPECT_FALSE(RoadTwinIndex::Load(bad, *reader, 30.0, 80.0, true, true, &why)) << "truncated";
  auto flipped = bytes;
  flipped[0] ^= 0x20;
  write_all(bad, flipped);
  EXPECT_FALSE(RoadTwinIndex::Load(bad, *reader, 30.0, 80.0, true, true, &why)) << "magic";

  // The engine's entry point scans instead of using a refused file.
  const auto& fallback = RoadTwinIndex::get(*other_reader, 30.0, 80.0, true, true, path);
  EXPECT_FALSE(fallback.mapped());
  EXPECT_TRUE(fallback.same_contents(*RoadTwinIndex::Scan(*other_reader, 30.0, 80.0, true, true)));
}

// ---------------------------------------------------------------------------------
// The engine, end to end, on a map whose loop depends on the index.  25 m/char.
//
//   A------------------------------------------------B     A-B curvature 15
//
//   P------------------------------------------------Q     50 m: the parallel street
//
//   F------------------------------------------------E     200 m: the fresh way home
//
// With the parallel tier the return leaves P-Q for E-F; without it, it rides P-Q
// (test_roundtrip_audit P1a). A mapped index that lost its parallels would serve the
// P-Q loop.
// ---------------------------------------------------------------------------------
class RoadIdentityPersistRoundTrip : public ::testing::Test {
protected:
  static gurka::map map;
  static void SetUpTestSuite() {
    const std::string ascii_map = R"(
    A----------------------------------------------------------------------------------------------------------------------B

    P----------------------------------------------------------------------------------------------------------------------Q




    F----------------------------------------------------------------------------------------------------------------------E
    )";
    const gurka::ways ways = {
        {"AB", {{"highway", "secondary"}}}, {"BQ", {{"highway", "secondary"}}},
        {"QP", {{"highway", "secondary"}}}, {"PA", {{"highway", "secondary"}}},
        {"BE", {{"highway", "secondary"}}}, {"EF", {{"highway", "secondary"}}},
        {"FA", {{"highway", "secondary"}}},
    };
    const auto layout = gurka::detail::map_to_coordinates(ascii_map, 25);
    map = gurka::buildtiles(layout, ways, {}, {}, "test/data/road_identity_persist_rt");
    map.config.put("thor.roundtrip_pair_pass", false);
    map.config.put("thor.roundtrip_xcand_penalty", false);
    map.config.put("thor.roundtrip_gate_refill_budget", 0);

    auto reader = test::make_clean_graphreader(map.config.get_child("mjolnir"));
    std::vector<baldr::GraphId> curvy;
    for (const auto& [a, b] :
         std::vector<std::pair<std::string, std::string>>{{"A", "B"}, {"B", "A"}})
      curvy.push_back(std::get<0>(gurka::findEdgeByNodes(*reader, layout, a, b)));
    test::customize_edges(map.config,
                          [&curvy](const baldr::GraphId& edgeid, baldr::DirectedEdge& edge) {
                            if (std::find(curvy.begin(), curvy.end(), edgeid) != curvy.end())
                              edge.set_curvature(15);
                          });
  }
  static std::string run(gurka::map m, const std::string& file, bool parallel_tier = true) {
    m.config.put("thor.roundtrip_parallel_tier", parallel_tier);
    m.config.put("thor.roundtrip_road_identity_file", file);
    std::string json;
    gurka::do_action(valhalla::Options::route, m, {"A", "A"}, "motorcycle",
                     {{"/roundtrip/target_distance", "6000"},
                      {"/roundtrip/num_candidates", "1"},
                      {"/costing_options/motorcycle/reuse_penalty", "0.8"},
                      {"/costing_options/motorcycle/prefer_curvature", "0.5"}},
                     {}, &json);
    return json;
  }
};
gurka::map RoadIdentityPersistRoundTrip::map = {};

TEST_F(RoadIdentityPersistRoundTrip, MappedIndexServesTheSameLoop) {
  const std::string path = file_beside(map, "road_identity.bin");
  {
    auto reader = test::make_clean_graphreader(map.config.get_child("mjolnir"));
    ASSERT_TRUE(RoadTwinIndex::Scan(*reader, 30.0, 80.0, true, true)->Save(path, *reader));
  }
  const std::string scanned = run(map, "");
  const std::string mapped = run(map, path);
  ASSERT_NE(scanned.find("\"trip\""), std::string::npos) << scanned;

  auto reader = test::make_clean_graphreader(map.config.get_child("mjolnir"));
  ASSERT_TRUE(RoadTwinIndex::get(*reader, 30.0, 80.0, true, true, path).mapped())
      << "the engine did not use the persisted file";
  EXPECT_EQ(scanned, mapped);

  // Control: the loop depends on the index, so equality above is not vacuous.
  EXPECT_NE(scanned, run(map, "", false)) << "map not discriminating";
}
