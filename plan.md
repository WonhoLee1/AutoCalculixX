
gmsh - 형상 및 메쉬 생성, 멀티 파트 가능, set(node, element)

calculix - 입력 파일 export (abaqus)

고유진동수 해석을 위한 calculix 입력 파일 수정

calulix 실행을 통한 결과 파일 (frd, dat) 생성

calulix 결과 파일 (dat) 분석을 통해 고유진동수 관련 정보 획득(출력)

frd to vtk 변환 스크립트 또는 명령 프로그램 실행

pyvista viewer를 간단하게 만들어서 vtk 보여주기 (mode shape)

체계적으로 확장 가능한 객체지향 구조 기반
임의의 inp (abaqus) 파일이 입력으로 들어오면, 고유진동수 해석을 할 수 있도록 하는 기능
inp파일이 include 되는 구조

파이프라인 테스트에는 다음의 소스 코드 참고, gmsh로 tray 형상을 만들고 mesh를 생성한 후 내부 솔버로 고유진동수를 수행
D:\PythonCodeStudy\WHT_LightChassisModel\test_jaxSSO\exam2_shell_jaxSSO.py

calculix 실행 파일 설치 위치
D:\SOFTWARE\calculix_2.23_4win

참고 사이트

<https://www.calculix.de/>

<https://github.com/mhayrettin/PrePoMax_CalculiX_OpenSource/blob/main/getResultsFromFRD/getResultsFromFRD.py>

<https://github.com/mhayrettin/PrePoMax_CalculiX_OpenSource/tree/main/getModalResultsFromDAT>

<https://github.com/rsmith-nl/calculix-frdconvert/blob/main/frdconvert.py>

<https://deepwiki.com/calculix/CalculiX-Examples/1.1-getting-started>

<https://deepwiki.com/calculix/CalculiX-Examples/3.7-modal-analysis>

먼저 구조부터 짜보자.
