{ fetchFromGitHub, python3Packages }:

python3Packages.buildPythonPackage {
  pname = "squmfit";
  version = "0.2";
  src = fetchFromGitHub {
    owner = "bgamari";
    repo = "squmfit";
    rev = "e1d14055eb8f840d837b7e10f11a8ddfdac208c7";
    sha256 = "0486rr4y8nik4mrpwlh2xswqr2gpiqwp7lpfj7jqjh68ycdlfk9h";
  };
  propagatedBuildInputs = with python3Packages; [numpy scipy matplotlib];
}

