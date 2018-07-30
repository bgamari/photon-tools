{ fetchFromGitHub, python3Packages }:

let
  squmfit = import ./squmfit.nix { inherit fetchFromGitHub python3Packages; };
in python3Packages.buildPythonPackage {
  pname = "photon-tools";
  version = "0.1";
  src = ./.;
  propagatedBuildInputs = with python3Packages; [numpy scipy matplotlib squmfit];
  nativeBuildInputs = with python3Packages; [cython];
}
