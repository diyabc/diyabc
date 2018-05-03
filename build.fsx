#r "paket:
nuget Fake.IO.FileSystem
nuget Fake.Core.Target //"
#load "./.fake/build.fsx/intellisense.fsx"

open Fake.Core
open Fake.IO
open Fake.IO.Globbing.Operators

let filterDir d f =
    f
    |> Seq.fold (fun s t -> s -- (d + t))
                !! (d + "*")

let testsnpdir = __SOURCE_DIRECTORY__ + @"/../F4_FD_Collin/f4mat_test_snp/"
let testsnpfiles = [ "header.txt"; "RNG_state_0000.bin"; "sim_dataset_001.snp" ]
Target.create "Test_snp" (fun _ ->
    (filterDir testsnpdir testsnpfiles)
    |> File.deleteAll
)

let testsnpheaderdir = __SOURCE_DIRECTORY__ + @"/../F4_FD_Collin/f4mat_header_snp/"
let testsnpheaderfiles = [ "header.txt"; "header-orig.txt"; "RNG_state_0000.bin"; "sim_dataset_001.snp" ]

Target.create "Header_snp" (fun _ ->
    (filterDir testsnpheaderdir testsnpheaderfiles)
    |> File.deleteAll
)

let testmatdir = __SOURCE_DIRECTORY__ + @"/../F4_FD_Collin/f4mat_test/"
let testmatfiles = [ "header.txt"; "RNG_state_0000.bin"; "exemple_dataset_poolseq.pool"; "xtest"]
Target.create "Test_mat" (fun _ ->
    (filterDir testmatdir testmatfiles)
    |> File.deleteAll
)

let testxdir = __SOURCE_DIRECTORY__ + @"/../F4_FD_Collin/f4mat_test/xtest/"

Target.create "Test_x" (fun _ ->
    !! (testxdir + "*")
    |> File.deleteAll
    Shell.cd (testxdir + "../")
    [ "exemple_dataset_poolseq.pool"
      "reftableRF.bin"
      "RNG_state_0000.bin"
      "headerRF.txt" ]
    |> Shell.copyFiles testxdir 
    Shell.cd testxdir
    Shell.rename "header.txt" "headerRF.txt"
    Shell.rename "reftable.bin" "reftableRF.bin"
    Shell.cd __SOURCE_DIRECTORY__
)

let rangertestbed = __SOURCE_DIRECTORY__ + @"../ranger-testbed/"

Target.create "Ranger_Testbed" (fun _ ->
    !! (rangertestbed + "ranger_*")
    |> File.deleteAll
)

Target.runOrDefault "Test_x"