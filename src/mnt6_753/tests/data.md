## for ec.rs

```
mnt6q = 418984909679189534023442147912406371281707099199539490717835029210253
....: 52812571106773058893763790338921418070971888458477323173057491593855069696241
....: 85479639616572141632535006444147041813784639846961193571905990816422078447616
....: 0001

mnt6a = 11
mnt6b = 116259089995413211520273402240103747168411677017835846483389082354108
....: 59267060079819722747939267925389062611062156601938166010098747920378738927832
....: 65813362545426011540907581618755505585949025337570472802794431550112272342687
....: 9114

mnt6 = EllipticCurve(GF(mnt6q), [mnt6a, mnt6b])
```


g1 is valid (just need to check is on curve)
```
u = (28090895838563800711133879843595433775537169677618939443143986918075792449676954374829149809740383781044105267184049529844255718616804686728808462794011444016618394095087707440046649986355174030866091888972335463663274809274579 : 4551579903480647055886250665303577983556258624650980342070285056070310929270852920443313052780891531023164520967524747091364975334493394664683782769215353073859327289253789705786348722457868585236746830711252356374189997310018 : 1)
```

```
u = (ux, uy)
ux = 0xa402e79179c2b0d3
0xd8ca333cff73ad6b
0x3576bcf14412822b
0x7e334a4abf49a2c
0xda4c0ce7eb3ec362
0xaab59aba7c0cf0ed
0x63aa8468ad19a19f
0xe4c86d0117c4cd21
0x96749f1c07e2fbac
0x57ecb3c5dcf96135
0x4b56cce8ffa4e36
0x12f901fd112e5
uy = 0x879d504f2cb84042
0x6e1466b9f8302dcc
0x28bb71efe4d3ec8a
0xeb78daf4ee8f0000
0x68e22597335f1092
0xb484c7db5dd1d830
0x11348d02eabbcd78
0x3181125f0a287b36
0xb706941e3cf9af1d
0xef396d80a3fb04a
0xf917d61a1f02eed3
0x312fba94f7ed

```



no subgroups

g1 addition
```
u = (ux : uy : 1) = (24606587867320141622595849793297499470128572816163062973210476902258668583052647838688576726142879981604823038579416607770453028821688898367467695712011719839151059828410346256782254827585447503860119492354187732503887980526633 : 21992562594513691142231533487896168949669612362414795017453481841069568342928897702139514232800985283895173271225210586933289850699985925577661765337058158738895419982898529913636922510358337210520991074924670517306030679832342 : 1)
v = (vx : vy : 1) = (39822029986824284207085257036873227523555264909356540756785454424475392883475582937924256490223668919767511726819573829332562386261930595064563237848509046876988509300761150674681298061692149487299749736708476502700329848655499 : 35219752738518272611000219044555180789086811598253750085684891442822511822472908789845263266400549421470137094556082171109352552066029863048127731701959909463757254336972300776886445168971710339109313957142267721327865337415057 : 1)
u + v = (wx : wy : 1) = (22461654423102492156408206221825670259975969584549351327649543791614885655594157022430358240129041696178453291328360413312604916409752339781109860389237710910044892697429514958071367569671470846588689572886546022131989630893701 : 23466963882044990799639566269946709879639178211180008608934074647442159332776972641137334405322403472972259780132498754972994981145273513090684325307727345871038982033708089606816662043822126227674564537690460454046917171721015 : 1)
```
where
```
ux = 0xfa3269471af7fc29
0xe75814faa6e78c3
0x49a7d9cd2dd6aa5b
0x78cb72dd40780e8f
0xa628467279c19e1
0xaa2dabe7aacacfa4
0xfd0449a399f5bc10
0xcc69871930e2b94
0xbd359add50190868
0xa430620337195647
0x870751af8ead4105
0x109e8f3cec114

uy = 0xc072a5045a0eb716
0xaa7ea85fd88b5c40
0x28ba48f5173291fa
0x15530eef5481573d
0xdd1ad67618327dac
0xf5c0666ed192a3b2
0x2b70dfe53a20a8ce
0x3d354613aff0bd6a
0x7a0034d2e30de04b
0xc3eabcb628fe0291
0x5b3196b3e0ce308a
0xeda96045ae66

vx = 0xb6f0a8e0b589ce8b
0x11d2007be6d43bcb
0x8cf5e84f9d1455e7
0xb3948e0450cfad43
0x3150075f6db70ae1
0x227254395f6e7f5c
0xd4013b9d00959599
0xb9b8b393ab316213
0x7cdd3853be3db699
0xb00bdce397d661ff
0xc1e075f0736f6f8d
0x1ae55bf692b13

vy = 0x39354a603130591
0x65ef5dc8a3bb6cad
0x3dc56a836204653f
0x189514b20bb6a1e
0x51d29d2c02d4e641
0x1f7a4d415b6be19a
0x144dc1b1ffcc41b0
0x485d891f0870162b
0xa8a4ca2c5c8d5dd0
0xbfb746c7be9ecbfa
0x2bd71fee34ba3197
0x17c99c45f0909

wx = 0x4075f3c009273e85
0x5d4bd293f6bcc9e6
0x8e27c491f7573c54
0xbd6387c831a6ca54
0xc2dca28819554d1e
0x347a4cfd65ecc90b
0xcf2e71ba8fa7bb4c
0x72338f027482c52c
0x870a69b5ee0b8e64
0xe54577a0ca314b35
0x4227f07bb68b1667
0xf2bb18a3958a

wy = 0x197f59872551b337
0xea7728119a1c4aee
0xa86d3ead69ab4e97
0x2310388bf7be5af3
0x3e2810110e19cfb2
0x71ff5fc55bf03cc9
0xe6b3608817306d89
0xeb20fd323630411d
0x7294ef5c1a0476b9
0x33ed661992c99ae2
0x593721f98fddb53d
0xfd983ca844c2
```


g1 doubling
```
u = (ux : uy : 1) = (31746117892822807594504339111303968665364263067619663603839573471912247332556604750632398214343784139578300305101183556892617284794055445777667973235741057454617109362922688706335579079869367210192811531861741378391352069414793 : 18558771584904833921658772569363477475155163851918659828453828825594397922771834506876543033540979054833419847438850768557405726732076226138088088364210265121673756914377049682210337038949635990064794019837554376934424045520875 : 1)
2u = (u2x : u2y : 1) = (27099887515698502538241628710808293449882448043133987854781216362006448082690394272453166898981529388881828135195114301143342488776643293517672429592805957233794254278966042853994035595216143754353538027408734210216612745320186 : 12471417379020085991740073782362257416401746533213966334096043008343713562421562315988728187554836788091965338672448203073765825388404631396125862255433594742883998897195193839102813036886630568540299626318478649109996362961104 : 1)
```
where
```
ux = 0xb19aed9b3dc15f89
0x29de7cb8bc671412
0x4c7542b9dea227d4
0x5dc4af9e491c7a5f
0xdee5857fa6bc0ca2
0xb26ed03d5845df9a
0x55dd26801730f0c7
0x277ca14f3ca740dc
0xcd92eab91d54a80b
0xb0511416f2a0a418
0x3002007605fcca03
0x157101f321264

uy = 0xf1b7d78c336187eb
0x759be3b9a52861ef
0xbf535e3141e47420
0xc88cb29d5479ce23
0x8048c51c4a374a81
0xad6e60b0a50a47a7
0x15705b02b85cad3f
0x477910253eb87de5
0xe8c2425b4c8d3b25
0x8e64e07d1d7ab34a
0x3a4865ebe3d0bc24
0xc88df4f8a73c

u2x = 0x532c93cfcb4a16fa
0x5114f616cce055a6
0xf314349b3be478bd
0x421f32dbdfefbd72
0xf77cce337c6eacad
0xee7b6035844afe4d
0x31c3609e58efa430
0xa135db8c9f8a40bf
0xa4e02fe151db05ee
0x6ab8eb33c99569d3
0xf610aaeea6d03ac5
0x124da8c081bc9

u2y = 0xf2e6d63d02d690d0
0x742166cca7aa2dd1
0xf878b8f47cd2b343
0xbb6c5aa896dba201
0x844652a50840d74c
0xfed24238534ec2be
0x32d3c6d354e0fec4
0xb30eb0ae74a3fe0a
0xa799f9a48370cded
0x7511b3ce4a912397
0x537b51139b80f575
0x86c594e940e2
```


g2 is valid (just need to check is on curve)
u =

reject point on twist?

no subgroups

g2 addition
u =
v =
u + v =

g2 doubling
u =
2u =