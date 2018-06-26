{ fetchFromGitHub, python3Packages }:

let
  squmfit = python3Packages.buildPythonPackage {
    pname = "squmfit";
    version = "0.1";
    src = fetchFromGitHub {
      owner = "bgamari";
      repo = "squmfit";
      rev = "e1d14055eb8f840d837b7e10f11a8ddfdac208c7";
      sha256 = "0486rr4y8nik4mrpwlh2xswqr2gpiqwp7lpfj7jqjh68ycdlfk9h";
    };
    propagatedBuildInputs = with python3Packages; [numpy scipy matplotlib];
  };
in python3Packages.buildPythonPackage {
  pname = "photon-tools";
  version = "0.1";
  src = ./.;
  propagatedBuildInputs = with python3Packages; [numpy scipy matplotlib squmfit];
  nativeBuildInputs = with python3Packages; [cython];
}
