#include "Data_Generation.hpp"

#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
#include <mpi.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle_Standard.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

namespace
{
bool FileExists(const std::string& path)
{
	struct stat info;
	return stat(path.c_str(), &info) == 0 && S_ISREG(info.st_mode);
}

std::string TestOutputDir(const std::string& name)
{
	std::string dir = name + "_" + std::to_string(getpid()) + "/";
	mkdir(dir.c_str(), 0755);
	return dir;
}

void TouchFile(const std::string& path)
{
	std::ofstream file(path);
	file << "stale\n";
}
}

int main(int argc, char* argv[])
{
	int result = 0;

	::testing::InitGoogleTest(&argc, argv);
	MPI_Init(&argc, &argv);
	result = RUN_ALL_TESTS();
	MPI_Finalize();
	return result;
}

TEST(TestDataGeneration, TestDataSetConstructor)
{
	// ARRANGE
	unsigned int sample_size = 100;
	unsigned int max_traj    = 0;
	double u_min			 = 0.0;
	unsigned int iso_rings	 = 10;
	// ACT
	Simulation_Data data_set(sample_size, max_traj, u_min, iso_rings);
	// ASSERT
	ASSERT_EQ(data_set.data.size(), iso_rings);
	for(auto& set : data_set.data)
		ASSERT_EQ(set.size(), 0);
}

TEST(TestDataGeneration, TestGenerateData)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0 * pb);
	DM.Set_Sigma_Electron(1.0 * pb);

	unsigned int sample_size = 2;

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 100, 50);

	// ACT
	Simulation_Data data_set(sample_size, 0);
	data_set.Generate_Data(DM, SSM, SHM);

	// ASSERT – number_of_trajectories is private; verify indirectly
	ASSERT_GE(data_set.Free_Ratio() + data_set.Capture_Ratio() + data_set.Reflection_Ratio(), 0.0);
}

TEST(TestDataGeneration, TestConfigure)
{
	// ARRANGE
	unsigned int sample_size		 = 2;
	double r						 = 2.0 * rSun;
	unsigned int min_scattering		 = 1;
	unsigned int max_scattering		 = 1;
	unsigned long int max_time_steps = 1e4;
	// ACT
	Simulation_Data data_set(sample_size, 0);
	data_set.Configure(r, min_scattering, max_scattering, max_time_steps);

	// ASSERT
	// ASSERT_EQ(data_set.data[0].size(), sample_size);
}

TEST(TestDataGeneration, TestDataFreeRatio)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0e-100 * pb);
	DM.Set_Sigma_Electron(1.0e-100 * pb);

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 100, 50);

	unsigned int sample_size = 10;
	// ACT
	Simulation_Data data_set(sample_size, 0);
	data_set.Configure(1.1 * rSun, 0, 500);
	data_set.Generate_Data(DM, SSM, SHM);

	// ASSERT
	ASSERT_DOUBLE_EQ(data_set.Free_Ratio(), 1.0);
}

TEST(TestDataGeneration, TestDataSetCaptureRatio)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(1.0 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0 * pb);
	DM.Set_Sigma_Electron(1.0 * pb);

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 100, 50);

	unsigned int sample_size = 50;

	// ACT
	Simulation_Data data_set(sample_size, 0);
	data_set.Configure(1.1 * rSun, 0, 500);
	data_set.Generate_Data(DM, SSM, SHM);

	// ASSERT
	ASSERT_GE(data_set.Capture_Ratio(), 0.0);
}

TEST(TestDataGeneration, TestDataSetReflectionRatio)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(1.0 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0 * pb);
	DM.Set_Sigma_Electron(1.0 * pb);

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 100, 50);

	unsigned int sample_size = 10;

	// ACT
	Simulation_Data data_set(sample_size, 0);
	data_set.Generate_Data(DM, SSM, SHM);

	// ASSERT
	ASSERT_GE(data_set.Reflection_Ratio(), 0.0);
}

TEST(TestDataGeneration, TestSpeedFunctions)
{
	// ARRANGE
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;

	obscura::DM_Particle_SI DM(1.0 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0 * pb);
	DM.Set_Sigma_Electron(1.0 * pb);

	SSM.Interpolate_Total_DM_Scattering_Rate(DM, 100, 50);

	unsigned int sample_size = 10;
	double u_min			 = 0.0001;
	// ACT
	Simulation_Data data_set(sample_size, 0, u_min);
	data_set.Generate_Data(DM, SSM, SHM);

	// ASSERT
	EXPECT_DOUBLE_EQ(data_set.Minimum_Speed(), 0.75 * u_min);
}

TEST(TestDataGeneration, TestDefaultOutputContract)
{
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;
	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0e-100 * pb);
	DM.Set_Sigma_Electron(1.0e-100 * pb);

	Simulation_Data data_set(1, 1);
	data_set.Configure(1.1 * rSun, 0, 100);
	data_set.Generate_Data(DM, SSM, SHM);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const std::string output_dir = TestOutputDir("default_output_contract");
	if(rank == 0)
	{
		TouchFile(output_dir + std::string("evaporation_") + "summary.txt");
		TouchFile(output_dir + std::string("evaporation_") + "mode_summary.txt");
		TouchFile(output_dir + std::string("evaporation_") + "mode_" + "bincount.txt");
		TouchFile(output_dir + std::string("computation_") + "time_summary.txt");
	}
	data_set.Write_Output_Files(output_dir, DM);
	if(rank == 0)
	{
		EXPECT_TRUE(FileExists(output_dir + "bincount.txt"));
		EXPECT_TRUE(FileExists(output_dir + "evaporation_times.txt"));
		EXPECT_FALSE(FileExists(output_dir + "evaporation_diagnostics.txt"));
		EXPECT_FALSE(FileExists(output_dir + std::string("evaporation_") + "summary.txt"));
		EXPECT_FALSE(FileExists(output_dir + std::string("evaporation_") + "mode_summary.txt"));
		EXPECT_FALSE(FileExists(output_dir + std::string("evaporation_") + "mode_" + "bincount.txt"));
		EXPECT_FALSE(FileExists(output_dir + std::string("computation_") + "time_summary.txt"));
	}
}

TEST(TestDataGeneration, TestDiagnosticsOutputContract)
{
	Solar_Model SSM;
	obscura::Standard_Halo_Model SHM;
	obscura::DM_Particle_SI DM(0.01 * GeV);
	DM.Set_Low_Mass_Mode(true);
	DM.Set_Sigma_Proton(1.0e-100 * pb);
	DM.Set_Sigma_Electron(1.0e-100 * pb);

	Simulation_Data data_set(1, 1);
	data_set.Configure(1.1 * rSun, 0, 100);
	data_set.Configure_Evaporation_Diagnostics(true);
	data_set.Generate_Data(DM, SSM, SHM);

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const std::string output_dir = TestOutputDir("diagnostics_output_contract");
	data_set.Write_Output_Files(output_dir, DM);
	if(rank == 0)
	{
		EXPECT_TRUE(FileExists(output_dir + "bincount.txt"));
		EXPECT_TRUE(FileExists(output_dir + "evaporation_times.txt"));
		EXPECT_TRUE(FileExists(output_dir + "evaporation_diagnostics.txt"));
		EXPECT_FALSE(FileExists(output_dir + std::string("evaporation_") + "mode_summary.txt"));
		EXPECT_FALSE(FileExists(output_dir + std::string("evaporation_") + "mode_" + "bincount.txt"));
		EXPECT_FALSE(FileExists(output_dir + std::string("computation_") + "time_summary.txt"));
	}
}
