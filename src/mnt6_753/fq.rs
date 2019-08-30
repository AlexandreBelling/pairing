use super::fq3::Fq3;
use ff::{Field, PrimeField, PrimeFieldRepr};

// A coefficient of MNT6-753 G1, 11.
pub const A_COEFF: Fq = Fq(FqRepr([
    5145524327033718740,
    14149824967095184544,
    5159730833497260295,
    3902941467692815387,
    15830098551216085679,
    8665641533746801158,
    17502192300007146323,
    14483698255198590748,
    546300946688995976,
    4331975528992054828,
    5311428878520309260,
    495362057711802,
]));

// B coefficient of MNT6-753 G1, 11625908999541321152027340224010374716841167701783584648338908235410859267060079819722747939267925389062611062156601938166010098747920378738927832658133625454260115409075816187555055859490253375704728027944315501122723426879114.
pub const B_COEFF: Fq = Fq(FqRepr([
    8828711393625909642,
    12722539140758597443,
    2303826860244282256,
    8063890988281098391,
    6269149169423748670,
    3425772737529456013,
    1457017085322601211,
    5177155908178255133,
    18057960053344868113,
    10481469207136524576,
    17888199912367160320,
    290288558853910,
]));

// Generator of G1
// X = 16364236387491689444759057944334173579070747473738339749093487337644739228935268157504218078126401066954815152892688541654726829424326599038522503517302466226143788988217410842672857564665527806044250003808514184274233938437290,
// Y = 4510127914410645922431074687553594593336087066778984214797709122300210966076979927285161950203037801392624582544098750667549188549761032654706830225743998064330900301346566408501390638273322467173741629353517809979540986561128,
// Z = 1.
pub const G1_GENERATOR_X: Fq = Fq(FqRepr([
    13679520475591086443,
    2885136257016027368,
    11066012770060586598,
    5703030954402099790,
    3190768979802679266,
    6582995058780951235,
    5128324295984695312,
    1016733259901042723,
    17321836123078146775,
    8899743831920183727,
    9524148587687792298,
    115184260559151,
]));

pub const G1_GENERATOR_Y: Fq = Fq(FqRepr([
    12838766614814920872,
    11227272155367303431,
    9667298518712466279,
    2515147449921847221,
    12901526649129916826,
    14043676281050491370,
    6526808858766215037,
    2720702495038125485,
    5084284651197010970,
    7160614230104399452,
    4659354110445863515,
    215718083283178,
]));

pub const ZERO: Fq = Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));
pub const ONE: Fq = Fq(FqRepr([
    13373016969058414402,
    5670427856875409064,
    11667651089292452217,
    1113053963617943770,
    12325313033510771412,
    11510260603202358114,
    3606323059104122008,
    6452324570546309730,
    4644558993695221281,
    1127165286758606988,
    10756108507984535957,
    135547536859714,
]));

// A coefficient of MNT6-753 G2 =
// mnt6753_twist_coeff_a = mnt6753_Fq3(mnt6753_Fq::zero(), mnt6753_Fq::zero(),
//                                  mnt6753_G1::coeff_a);
//  = (ZERO, ZERO, A_COEFF);
pub const G2_A_COEFF: Fq3 = Fq3 {
    c0: ZERO,
    c1: ZERO,
    c2: A_COEFF,
};
// B coefficient of MNT6-753 G2 =
// mnt6753_twist_coeff_b = mnt6753_Fq3(mnt6753_G1::coeff_b * mnt6753_Fq3::non_residue,
//                                  mnt6753_Fq::zero(), mnt6753_Fq::zero());
// non_residue = mnt6753_Fq3::non_residue = mnt6753_Fq("11");
//  = (G1_B_COEFF * NON_RESIDUE, ZERO, ZERO);
//  =
//  (2189526091197672465268098090392210500740714959757583916377481826443393499947557697773546040576162515434508768057245887856591913752342600919117433675080691499697020523783784738694360040853591723916201150207746019687604267190251,
//  0, 0)
pub const G2_B_COEFF: Fq3 = Fq3 {
    c0: Fq(FqRepr([
        3284231658830416104,
        13720030246451177991,
        6276939417009443243,
        8340612253649729185,
        4863511590806861670,
        15883218135158530927,
        4865336109262680856,
        16600307443495218926,
        10112528487499131659,
        17308657107605697754,
        5326857497786417651,
        206191604157846,
    ])),
    c1: ZERO,
    c2: ZERO,
};


// From Coda/libff: non_residue = 11
pub const NON_RESIDUE: Fq = Fq(FqRepr([
    5145524327033718740,
    14149824967095184544,
    5159730833497260295,
    3902941467692815387,
    15830098551216085679,
    8665641533746801158,
    17502192300007146323,
    14483698255198590748,
    546300946688995976,
    4331975528992054828,
    5311428878520309260,
    495362057711802,
]));

// Generator of G2
// These are three Fq elements each because X and Y (and Z) are elements of Fq^3
// X = 46538297238006280434045879335349383221210789488441126073640895239023832290080310125413049878152095926176013036314720850781686614265244307536450228450615346834324267478485994670716807428718518299710702671895190475661871557310,
// 10329739935427016564561842963551883445915701424214177782911128765230271790215029185795830999583638744119368571742929964793955375930677178544873424392910884024986348059137449389533744851691082159233065444766899262771358355816328,
// 19962817058174334691864015232062671736353756221485896034072814261894530786568591431279230352444205682361463997175937973249929732063490256813101714586199642571344378012210374327764059557816647980334733538226843692316285591005879,
// Y = 5648166377754359996653513138027891970842739892107427747585228022871109585680076240624013411622970109911154113378703562803827053335040877618934773712021441101121297691389632155906182656254145368668854360318258860716497525179898,
// 26817850356025045630477313828875808893994935265863280918207940412617168254772789578700316551065949899971937475487458539503514034928974530432009759562975983077355912050606509147904958229398389093697494174311832813615564256810453,
// 32332319709358578441696731586704495581796858962594701633932927358040566210788542624963749336109940335257143899293177116050031684054348958813290781394131284657165540476824211295508498842102093219808642563477603392470909217611033,
// Z = 1.

pub const G2_GENERATOR_X_C0: Fq = Fq(FqRepr([
    10851632081534502623,
    2293958624975688101,
    188919536117884749,
    17860331619214007027,
    1887803711364835009,
    12024909542936874634,
    10408468762840763117,
    2218994925646485884,
    5852783618064312507,
    5105630124344468993,
    15624271947652004457,
    343291609063118,
]));

pub const G2_GENERATOR_X_C1: Fq = Fq(FqRepr([
    5125421133817098708,
    4635671833371751021,
    5021088233911155688,
    1910278417010696574,
    15893067539329770376,
    13275414949851777243,
    5574902941155169916,
    11685892668191778913,
    39603799641674991,
    15923769296743881818,
    2229648642576712242,
    271230940026616,
]));

pub const G2_GENERATOR_X_C2: Fq = Fq(FqRepr([
    13312714422637623631,
    7908936083790331333,
    16404290491319144526,
    10468940401598991065,
    16719607066400875106,
    12786278203212966733,
    12073070934010646425,
    3586662858255838647,
    15314405800929445490,
    17316524929057257003,
    11733547285987599962,
    369098436863856,
]));

pub const G2_GENERATOR_Y_C0: Fq = Fq(FqRepr([
    7838343860819620738,
    6475978215607342230,
    6907947951617093760,
    4877650884610899013,
    5572804598516141830,
    2857981775984637133,
    7047604422681304779,
    14272181022877590717,
    11871160228234478242,
    6227559544769149955,
    13648209694541221286,
    92423087846120,
]));

pub const G2_GENERATOR_Y_C1: Fq = Fq(FqRepr([
    7399349695289784545,
    10343186604152989839,
    5569933740221267149,
    16274755507370250244,
    17881568224199186382,
    13445460868909972033,
    8855491571276536451,
    4029003301327442037,
    17200334844901331034,
    5623722999341396242,
    1043926657758222554,
    345951167723219,
]));

pub const G2_GENERATOR_Y_C2: Fq = Fq(FqRepr([
    16283109937764686785,
    12936144524232368259,
    851410484573527812,
    13652212824151377038,
    1687600484300677308,
    13313718100971044033,
    7296597353247060371,
    11678393821640342983,
    3610823548198042284,
    9307794189036495061,
    7861574050710227414,
    380700544770643,
]));

// Coefficients for the Frobenius automorphism.
// c1[0] = 1,
// c1[1] = 24129022407817241407134263419936114379815707076943508280977368156625538709102831814843582780138963119807143081677569721953561801075623741378629346409604471234573396989178424163772589090105392407118197799904755622897541183052132
// c1[2] = 17769468560101711995209951371304522748355002843010440790806134764399814103468274958215310983651375801610927890210888755369611256415970113691066895445191924931148019336171640277697829047741006062493737919155152541323243293107868,
// c2 = {c1[0], c1[2], c1[1]}

// TODO: check if the values are correct
pub const FROBENIUS_COEFF_FQ3_C1: [Fq; 3] = [
    ONE,
    Fq(FqRepr([
        7739145380395648640,
        1403348385939055902,
        11220424057264707228,
        4567962295300549271,
        5929583493640677751,
        17618207486530478833,
        16600462137977359741,
        16551719371247820635,
        12057922785354578416,
        13022559182829558162,
        13308285686168533250,
        313705269181021,
    ])),
    Fq(FqRepr([
        12973180669431253567,
        17038664486452692616,
        11034024317238370177,
        7712681843988565810,
        4725787734130647531,
        2175028350442404679,
        9323639551697167751,
        14465264105466053583,
        8569442212929419360,
        17553812953652473294,
        13991744086792172309,
        48577617831792,
    ])),
];

// TODO: find value
// mnt6753_Fq3::Frobenius_coeffs_c2[0] = mnt6753_Fq("1");
// mnt6753_Fq3::Frobenius_coeffs_c2[1] = mnt6753_Fq("17769468560101711995209951371304522748355002843010440790806134764399814103468274958215310983651375801610927890210888755369611256415970113691066895445191924931148019336171640277697829047741006062493737919155152541323243293107868");
// mnt6753_Fq3::Frobenius_coeffs_c2[2] = mnt6753_Fq(24129022407817241407134263419936114379815707076943508280977368156625538709102831814843582780138963119807143081677569721953561801075623741378629346409604471234573396989178424163772589090105392407118197799904755622897541183052132");

pub const FROBENIUS_COEFF_FQ3_C2: [Fq; 3] = [
    ONE,
    Fq(FqRepr([
        12973180669431253567,
        17038664486452692616,
        11034024317238370177,
        7712681843988565810,
        4725787734130647531,
        2175028350442404679,
        9323639551697167751,
        14465264105466053583,
        8569442212929419360,
        17553812953652473294,
        13991744086792172309,
        48577617831792,
    ])),
    Fq(FqRepr([
        7739145380395648640,
        1403348385939055902,
        11220424057264707228,
        4567962295300549271,
        5929583493640677751,
        17618207486530478833,
        16600462137977359741,
        16551719371247820635,
        12057922785354578416,
        13022559182829558162,
        13308285686168533250,
        313705269181021,
    ])),
];

// Coefficients for the Frobenius automorphism.
// c1[0] = 1,
// c1[1] = 24129022407817241407134263419936114379815707076943508280977368156625538709102831814843582780138963119807143081677569721953561801075623741378629346409604471234573396989178424163772589090105392407118197799904755622897541183052133
// c1[2] = 24129022407817241407134263419936114379815707076943508280977368156625538709102831814843582780138963119807143081677569721953561801075623741378629346409604471234573396989178424163772589090105392407118197799904755622897541183052132
// c1[3] = 41898490967918953402344214791240637128170709919953949071783502921025352812571106773058893763790338921418070971888458477323173057491593855069696241854796396165721416325350064441470418137846398469611935719059908164220784476160000
// c1[4] = 17769468560101711995209951371304522748355002843010440790806134764399814103468274958215310983651375801610927890210888755369611256415970113691066895445191924931148019336171640277697829047741006062493737919155152541323243293107868
// c1[5] = 17769468560101711995209951371304522748355002843010440790806134764399814103468274958215310983651375801610927890210888755369611256415970113691066895445191924931148019336171640277697829047741006062493737919155152541323243293107869
pub const FROBENIUS_COEFF_FQ6_C1: [Fq; 6] = [
    ONE,
    Fq(FqRepr([
        2665418275744511426,
        7073776242814464967,
        4441331072847607829,
        5681016258918493042,
        18254896527151449163,
        10681724016023285331,
        1760041123371930134,
        4557299868084578750,
        16702481779049799698,
        14149724469588165150,
        5617650120443517591,
        449252806040736,
    ])),
    Fq(FqRepr([
        7739145380395648640,
        1403348385939055902,
        11220424057264707228,
        4567962295300549271,
        5929583493640677751,
        17618207486530478833,
        16600462137977359741,
        16551719371247820635,
        12057922785354578416,
        13022559182829558162,
        13308285686168533250,
        313705269181021,
    ])),
    Fq(FqRepr([
        2265581976117350591,
        18442012872391748519,
        3807704300793525789,
        12280644139289115082,
        10655371227771325282,
        1346491763263331896,
        7477357615964975877,
        12570239403004322603,
        2180620924574446161,
        12129628062772479841,
        8853285699251153944,
        362282887012814,
    ])),
    Fq(FqRepr([
        12973180669431253567,
        17038664486452692616,
        11034024317238370177,
        7712681843988565810,
        4725787734130647531,
        2175028350442404679,
        9323639551697167751,
        14465264105466053583,
        8569442212929419360,
        17553812953652473294,
        13991744086792172309,
        48577617831792,
    ])),
    Fq(FqRepr([
        7899453564780116353,
        4262348269618550065,
        4254931332821270779,
        8825735807606509581,
        17051100767641418943,
        13685288953644762793,
        12929962610801289759,
        2470844602302811697,
        13214001206624640642,
        234234166701528666,
        6301108521067156651,
        184125154691507,
    ])),
];

#[allow(unused)]
pub const NEGATIVE_ONE: Fq = Fq(FqRepr([
    2265581976117350591,
    18442012872391748519,
    3807704300793525789,
    12280644139289115082,
    10655371227771325282,
    1346491763263331896,
    7477357615964975877,
    12570239403004322603,
    2180620924574446161,
    12129628062772479841,
    8853285699251153944,
    362282887012814,
]));

/* Fq3 constants */

pub const FQ3_NQR_T: Fq3 = Fq3 {
    c0: Fq(FqRepr([
        2456656400918202012,
        7503386575313625620,
        1014314685003569848,
        10473903647598823719,
        15893393002146336511,
        8418203974290622500,
        9017296731996077946,
        2923126592994124774,
        9368756030960215800,
        17344552888362241070,
        10938255746876359306,
        107029542386399,
    ])),
    c1: ZERO,
    c2: ZERO,
};

pub const FQ3_T_MINUS_1: [u64; 36] = [
    15439605736802142541,
    18190868848461853149,
    6220121510046940818,
    10310485528612680366,
    5032137869959796540,
    3943048799800510054,
    1971151279016362045,
    6096644900171872841,
    12908407994230849218,
    4163225373804228290,
    10382959950522770522,
    9008828410264446883,
    18411821899404157689,
    12386199240837247984,
    13370099281150720481,
    11909278545073807560,
    5964354403900302648,
    15347506722065009035,
    7045354120681109597,
    14294096902719509929,
    6180325033003959541,
    14381489272445870003,
    18159920240207503954,
    17487026929061632528,
    12314108197538755669,
    12116872703077811769,
    3401400733784294722,
    13905351619889935522,
    10972472942574358218,
    6104159581753028261,
    4690139121547787552,
    4880965491878697414,
    1926648890365125214,
    13532564555356297305,
    3114545746551080,
    0,
];
/* instantiate */

#[derive(PrimeField)]
#[PrimeFieldModulus = "41898490967918953402344214791240637128170709919953949071783502921025352812571106773058893763790338921418070971888458477323173057491593855069696241854796396165721416325350064441470418137846398469611935719059908164220784476160001"]
#[PrimeFieldGenerator = "17"]
pub struct Fq(FqRepr);

impl Fq {
    pub fn mul_by_nonresidue(&mut self) {
        self.mul_assign(&NON_RESIDUE)
    }
}

pub const TWIST: Fq3 = Fq3 {
    c0: ZERO,
    c1: ONE,
    c2: ZERO,
};

pub const TWIST_INV: Fq3 = Fq3 {
    c0: ZERO,
    c1: ZERO,
    c2: Fq(FqRepr([
        7047753277478936932, 
        2312173296821124137, 
        10098270251729416131, 
        10600299174452634607, 
        1888293341983374451, 
        9380390556403650454, 
        15434851093720341428, 
        16202838202940280131, 
        11657871854320050050, 
        5446192953384660801, 
        3711759446995951857, 
        464895615962273,
    ])),
};

// 204691208819330962009469868104636132783269696790011977400223898462431810102935615891307667367766898917669754470400
pub const MNT6_X: [u64; 6] = [
    0x7a7713041ba18000,
    0x6b0344c4e2c428b0,
    0x733b714aa43c31a6,
    0x51852c8cbe26e600,
    0x86dcbcee5dcda7fe,
    0x015474b1d641a3fd,
];

pub const MNT6_X_IS_NEGATIVE: bool = false;

// 204691208819330962009469868104636132783269696790011977400223898462431810102935615891307667367766898917669754470400
pub const EXP_W0: [u64; 6] = [
    0x7a7713041ba18000,
    0x6b0344c4e2c428b0,
    0x733b714aa43c31a6,
    0x51852c8cbe26e600,
    0x86dcbcee5dcda7fe,
    0x015474b1d641a3fd,
];

pub const EXP_W0_IS_NEGATIVE: bool = false;

pub const EXP_W1: [u64; 1] = [1u64];

#[test]
fn test_a_coeff() {
    assert_eq!(Fq::from_repr(FqRepr::from(11)).unwrap(), A_COEFF);
}

#[allow(dead_code)]
fn print_repr(name: &str, s: &str) {
    let x = Fq::from_str(s).unwrap();
    match x {
        Fq(FqRepr(a)) => println!(
            "pub const {} : Fq = Fq(FqRepr([{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {},]));",
            name, a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11]
        ),
    }
}

#[test]
fn test_b_coeff() {
    
    // print_repr("a", "24670256968250050314369453338414040398886372940498750317546990920176458455051782638050675490736153330004092713609661781126033031593433591653542680180442304357679976011724949176664866695916116063153824205881192201645178417332921");
    // print_repr("b", "156791020203366180260256657231874686737199651945493055636390214454188587104614874034858177945948952768390490226316082318278825142036422930729931086032350144661055785447111343487177422537149755315674465512716951020508718684077");
    // print_repr("c", "15628311079616238046013893117672780759669673818011272427655444626315957136145939417351843067825738872121349656754567187264869708448001629773231887546884499289322960011567397346119732459010602759015282847807274942533131290273441");
    // print_repr("A_COEFF", "11");
    // print_repr("B_COEFF", "11625908999541321152027340224010374716841167701783584648338908235410859267060079819722747939267925389062611062156601938166010098747920378738927832658133625454260115409075816187555055859490253375704728027944315501122723426879114");
    // print_repr("G1_GENERATOR_X", "16364236387491689444759057944334173579070747473738339749093487337644739228935268157504218078126401066954815152892688541654726829424326599038522503517302466226143788988217410842672857564665527806044250003808514184274233938437290");
    // print_repr("G1_GENERATOR_Y", "4510127914410645922431074687553594593336087066778984214797709122300210966076979927285161950203037801392624582544098750667549188549761032654706830225743998064330900301346566408501390638273322467173741629353517809979540986561128");
    // print_repr("ZERO", "0");
    // print_repr("ONE", "1");
    // print_repr("G2_B_COEFF_C0", "2189526091197672465268098090392210500740714959757583916377481826443393499947557697773546040576162515434508768057245887856591913752342600919117433675080691499697020523783784738694360040853591723916201150207746019687604267190251");
    // print_repr("G2_GENERATOR_X_C0", "46538297238006280434045879335349383221210789488441126073640895239023832290080310125413049878152095926176013036314720850781686614265244307536450228450615346834324267478485994670716807428718518299710702671895190475661871557310");
    // print_repr("G2_GENERATOR_X_C1", "10329739935427016564561842963551883445915701424214177782911128765230271790215029185795830999583638744119368571742929964793955375930677178544873424392910884024986348059137449389533744851691082159233065444766899262771358355816328");
    // print_repr("G2_GENERATOR_X_C2", "19962817058174334691864015232062671736353756221485896034072814261894530786568591431279230352444205682361463997175937973249929732063490256813101714586199642571344378012210374327764059557816647980334733538226843692316285591005879");
    // print_repr("G2_GENERATOR_Y_C0", "5648166377754359996653513138027891970842739892107427747585228022871109585680076240624013411622970109911154113378703562803827053335040877618934773712021441101121297691389632155906182656254145368668854360318258860716497525179898");
    // print_repr("G2_GENERATOR_Y_C1", "26817850356025045630477313828875808893994935265863280918207940412617168254772789578700316551065949899971937475487458539503514034928974530432009759562975983077355912050606509147904958229398389093697494174311832813615564256810453");
    // print_repr("G2_GENERATOR_Y_C2", "32332319709358578441696731586704495581796858962594701633932927358040566210788542624963749336109940335257143899293177116050031684054348958813290781394131284657165540476824211295508498842102093219808642563477603392470909217611033");

    // print_repr("FROBENIUS_FQ3_COEFF_1", "24129022407817241407134263419936114379815707076943508280977368156625538709102831814843582780138963119807143081677569721953561801075623741378629346409604471234573396989178424163772589090105392407118197799904755622897541183052132");
    // print_repr("FROBENIUS_FQ3_COEFF_2", "17769468560101711995209951371304522748355002843010440790806134764399814103468274958215310983651375801610927890210888755369611256415970113691066895445191924931148019336171640277697829047741006062493737919155152541323243293107868");

    // print_repr("FROBENIUS_COEFF_FQ6_C1_1", "24129022407817241407134263419936114379815707076943508280977368156625538709102831814843582780138963119807143081677569721953561801075623741378629346409604471234573396989178424163772589090105392407118197799904755622897541183052133");
    // print_repr("FROBENIUS_COEFF_FQ6_C1_2", "24129022407817241407134263419936114379815707076943508280977368156625538709102831814843582780138963119807143081677569721953561801075623741378629346409604471234573396989178424163772589090105392407118197799904755622897541183052132");
    // print_repr("FROBENIUS_COEFF_FQ6_C1_3", "41898490967918953402344214791240637128170709919953949071783502921025352812571106773058893763790338921418070971888458477323173057491593855069696241854796396165721416325350064441470418137846398469611935719059908164220784476160000");
    // print_repr("FROBENIUS_COEFF_FQ6_C1_4", "17769468560101711995209951371304522748355002843010440790806134764399814103468274958215310983651375801610927890210888755369611256415970113691066895445191924931148019336171640277697829047741006062493737919155152541323243293107868");
    // print_repr("FROBENIUS_COEFF_FQ6_C1_5", "17769468560101711995209951371304522748355002843010440790806134764399814103468274958215310983651375801610927890210888755369611256415970113691066895445191924931148019336171640277697829047741006062493737919155152541323243293107869");

    print_repr("TWIST_INV", "19044768621781342455611006723291198694623049963615431396265228600466069460259593987754042619904699509735486805403844762419624117041633570486225564479452907348055189238795483837032008244475635668005425326845412801918538398254546");

    // print_repr("NEGATIVE_ONE", "41898490967918953402344214791240637128170709919953949071783502921025352812571106773058893763790338921418070971888458477323173057491593855069696241854796396165721416325350064441470418137846398469611935719059908164220784476160000");

    // print_repr("q-3/4", "31423868225939215051758161093430477846128032439965461803837627190769014609428330079794170322842754191063553228916343857992379793118695391302272181391097297124291062244012548331102813603384798852208951789294931123165588357120000");
    // print_repr("q-1/2", "20949245483959476701172107395620318564085354959976974535891751460512676406285553386529446881895169460709035485944229238661586528745796927534848120927398198082860708162675032220735209068923199234805967859529954082110392238080000");
    
    // print_repr("root of unity",
    //     "5431548564651772770863376209190533321743766006080874345421017090576169920304713950094628043692772801995471539849411522704471393987882883355624697206026582300050878644000631322086989454860102191886653186986980927065212650747291");
    // print_repr("nqr", "22168644070733283197994897338612733221095941481265408161807376791727499343083607817089033595478370212662133368413166734396127674284827734481031659015434501966360165723728649019457855887066657739809176476252080335185730833468062");
    // assert_eq!(Fq::from_str("11625908999541321152027340224010374716841167701783584648338908235410859267060079819722747939267925389062611062156601938166010098747920378738927832658133625454260115409075816187555055859490253375704728027944315501122723426879114").unwrap(), B_COEFF);
}

#[cfg(test)]
use rand::{Rand, SeedableRng, XorShiftRng};

#[test]
fn test_neg_one() {
    let mut o = Fq::one();
    o.negate();

    assert_eq!(NEGATIVE_ONE, o);
}

#[test]
fn test_fq_repr_from() {
    assert_eq!(
        FqRepr::from(100),
        FqRepr([100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    );
}

#[test]
fn test_fq_repr_is_odd() {
    assert!(!FqRepr::from(0).is_odd());
    assert!(FqRepr::from(0).is_even());
    assert!(FqRepr::from(1).is_odd());
    assert!(!FqRepr::from(1).is_even());
    assert!(!FqRepr::from(324834872).is_odd());
    assert!(FqRepr::from(324834872).is_even());
    assert!(FqRepr::from(324834873).is_odd());
    assert!(!FqRepr::from(324834873).is_even());
}

#[test]
fn test_fq_repr_is_zero() {
    assert!(FqRepr::from(0).is_zero());
    assert!(!FqRepr::from(1).is_zero());
    assert!(!FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]).is_zero());
}

#[test]
fn test_fq_repr_sub_noborrow() {
    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    for _ in 0..1000 {
        let mut a = FqRepr::rand(&mut rng);
        a.0[11] >>= 30;
        let mut b = a;
        for _ in 0..10 {
            b.mul2();
        }
        let mut c = b;
        for _ in 0..10 {
            c.mul2();
        }

        assert!(a < b);
        assert!(b < c);

        let mut csub_ba = c;
        csub_ba.sub_noborrow(&b);
        csub_ba.sub_noborrow(&a);

        let mut csub_ab = c;
        csub_ab.sub_noborrow(&a);
        csub_ab.sub_noborrow(&b);

        assert_eq!(csub_ab, csub_ba);
    }
}

#[test]
fn test_fq_repr_add_nocarry() {
    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    // Test for the associativity of addition.
    for _ in 0..1000 {
        let mut a = FqRepr::rand(&mut rng);
        let mut b = FqRepr::rand(&mut rng);
        let mut c = FqRepr::rand(&mut rng);

        // Unset the first few bits, so that overflow won't occur.
        a.0[5] >>= 3;
        b.0[5] >>= 3;
        c.0[5] >>= 3;

        let mut abc = a;
        abc.add_nocarry(&b);
        abc.add_nocarry(&c);

        let mut acb = a;
        acb.add_nocarry(&c);
        acb.add_nocarry(&b);

        let mut bac = b;
        bac.add_nocarry(&a);
        bac.add_nocarry(&c);

        let mut bca = b;
        bca.add_nocarry(&c);
        bca.add_nocarry(&a);

        let mut cab = c;
        cab.add_nocarry(&a);
        cab.add_nocarry(&b);

        let mut cba = c;
        cba.add_nocarry(&b);
        cba.add_nocarry(&a);

        assert_eq!(abc, acb);
        assert_eq!(abc, bac);
        assert_eq!(abc, bca);
        assert_eq!(abc, cab);
        assert_eq!(abc, cba);
    }
}

#[test]
fn test_fq_is_valid() {
    let mut a = Fq(MODULUS);
    assert!(!a.is_valid());
    a.0.sub_noborrow(&FqRepr::from(1));
    assert!(a.is_valid());
    assert!(Fq(FqRepr::from(0)).is_valid());

    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    for _ in 0..1000 {
        let a = Fq::rand(&mut rng);
        assert!(a.is_valid());
    }
}

#[test]
fn test_fq_add_assign() {
    {
        // Random number
        let mut tmp = Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));
        assert!(tmp.is_valid());
        // Test that adding zero has no effect.
        tmp.add_assign(&Fq(FqRepr::from(0)));
        assert_eq!(tmp, Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,])));
        // Add one and test for the result.
        tmp.add_assign(&Fq(FqRepr::from(0)));
        assert_eq!(tmp, Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,])));
        // Add another random number that exercises the reduction.
        tmp.add_assign(&Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])));
        assert_eq!(tmp, Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,])));
        // Add one to (q - 1) and test for the result.
        tmp = Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));
        tmp.add_assign(&Fq(FqRepr::from(0)));
        assert!(tmp.0.is_zero());
        // Add a random number to another one such that the result is q - 1
        tmp = Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));
        tmp.add_assign(&Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])));
        assert_eq!(tmp, Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,])));
        // Add one to the result and test for it.
        tmp.add_assign(&Fq(FqRepr::from(0)));
        assert!(tmp.0.is_zero());
    }
    // Test associativity

    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    for _ in 0..1000 {
        // Generate a, b, c and ensure (a + b) + c == a + (b + c).
        let a = Fq::rand(&mut rng);
        let b = Fq::rand(&mut rng);
        let c = Fq::rand(&mut rng);

        let mut tmp1 = a;
        tmp1.add_assign(&b);
        tmp1.add_assign(&c);

        let mut tmp2 = b;
        tmp2.add_assign(&c);
        tmp2.add_assign(&a);

        assert!(tmp1.is_valid());
        assert!(tmp2.is_valid());
        assert_eq!(tmp1, tmp2);
    }
}

#[test]
fn test_fq_sub_assign() {
    {
        // Test arbitrary subtraction that tests reduction.
        let mut tmp = Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));
        tmp.sub_assign(&Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])));
        assert_eq!(tmp, Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,])));

        // Test the opposite subtraction which doesn't test reduction.
        tmp = Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));
        tmp.sub_assign(&Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])));
        assert_eq!(tmp, Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,])));

        // Test for sensible results with zero
        tmp = Fq(FqRepr::from(0));
        tmp.sub_assign(&Fq(FqRepr::from(0)));
        assert!(tmp.is_zero());

        tmp = Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));
        tmp.sub_assign(&Fq(FqRepr::from(0)));
        assert_eq!(tmp, Fq(FqRepr([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,])));
    }
    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    for _ in 0..1000 {
        // Ensure that (a - b) + (b - a) = 0.
        let a = Fq::rand(&mut rng);
        let b = Fq::rand(&mut rng);

        let mut tmp1 = a;
        tmp1.sub_assign(&b);

        let mut tmp2 = b;
        tmp2.sub_assign(&a);

        tmp1.add_assign(&tmp2);
        assert!(tmp1.is_zero());
    }
}

#[test]
fn test_fq_mul_assign() {
    let mut tmp = Fq(FqRepr([
        5910615766705689880,
        5941843556834618629,
        17738392444264801100,
        8519155207957497272,
        5644166707453047522,
        6508585849126723056,
        11551451519750857855,
        4923836763032933374,
        18200608768347826728,
        17755494014637076452,
        14535888706801877344,
        404104852082166,
    ]));
    tmp.mul_assign(&Fq(FqRepr([
        1491882177966700990,
        13517713027520413973,
        1434771261063841555,
        13808524160895852674,
        17216040672255697988,
        16728021351469661868,
        12143441390778418870,
        15190938871306883424,
        16465767661012805748,
        4524186958912820603,
        6884901866695231333,
        474902058234092,
    ])));
    assert!(
        tmp == Fq(FqRepr([
            11278943521444194658,
            1314698642669499282,
            4149508191745079138,
            5174227224058824588,
            17318046714871393845,
            976602439406741168,
            55276789091328686,
            1921588846044492479,
            15523310107835249280,
            5955525442673681304,
            10490137305032061504,
            327249176020484,
        ]))
    );

    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    for _ in 0..1000000 {
        // Ensure that (a * b) * c = a * (b * c)
        let a = Fq::rand(&mut rng);
        let b = Fq::rand(&mut rng);
        let c = Fq::rand(&mut rng);

        let mut tmp1 = a;
        tmp1.mul_assign(&b);
        tmp1.mul_assign(&c);

        let mut tmp2 = b;
        tmp2.mul_assign(&c);
        tmp2.mul_assign(&a);

        assert_eq!(tmp1, tmp2);
    }

    for _ in 0..1000000 {
        // Ensure that r * (a + b + c) = r*a + r*b + r*c

        let r = Fq::rand(&mut rng);
        let mut a = Fq::rand(&mut rng);
        let mut b = Fq::rand(&mut rng);
        let mut c = Fq::rand(&mut rng);

        let mut tmp1 = a;
        tmp1.add_assign(&b);
        tmp1.add_assign(&c);
        tmp1.mul_assign(&r);

        a.mul_assign(&r);
        b.mul_assign(&r);
        c.mul_assign(&r);

        a.add_assign(&b);
        a.add_assign(&c);

        assert_eq!(tmp1, a);
    }
}

#[test]
fn test_fq_squaring() {
    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    for _ in 0..1000000 {
        // Ensure that (a * a) = a^2
        let a = Fq::rand(&mut rng);

        let mut tmp = a;
        tmp.square();

        let mut tmp2 = a;
        tmp2.mul_assign(&a);

        assert_eq!(tmp, tmp2);
    }
}

#[test]
fn test_fq_inverse() {
    assert!(Fq::zero().inverse().is_none());

    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    let one = Fq::one();

    for _ in 0..1000 {
        // Ensure that a * a^-1 = 1
        let mut a = Fq::rand(&mut rng);
        let ainv = a.inverse().unwrap();
        a.mul_assign(&ainv);
        assert_eq!(a, one);
    }
}

#[test]
fn test_fq_double() {
    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    for _ in 0..1000 {
        // Ensure doubling a is equivalent to adding a to itself.
        let mut a = Fq::rand(&mut rng);
        let mut b = a;
        b.add_assign(&a);
        a.double();
        assert_eq!(a, b);
    }
}

#[test]
fn test_fq_negate() {
    {
        let mut a = Fq::zero();
        a.negate();

        assert!(a.is_zero());
    }

    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    for _ in 0..1000 {
        // Ensure (a - (-a)) = 0.
        let mut a = Fq::rand(&mut rng);
        let mut b = a;
        b.negate();
        a.add_assign(&b);

        assert!(a.is_zero());
    }
}

#[test]
fn test_fq_pow() {
    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    for i in 0..1000 {
        // Exponentiate by various small numbers and ensure it consists with repeated
        // multiplication.
        let a = Fq::rand(&mut rng);
        let target = a.pow(&[i]);
        let mut c = Fq::one();
        for _ in 0..i {
            c.mul_assign(&a);
        }
        assert_eq!(c, target);
    }

    for _ in 0..1000 {
        // Exponentiating by the modulus should have no effect in a prime field.
        let a = Fq::rand(&mut rng);

        assert_eq!(a, a.pow(Fq::char()));
    }
}

#[test]
fn test_fq_sqrt() {
    use ff::SqrtField;

    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    assert_eq!(Fq::zero().sqrt().unwrap(), Fq::zero());

    for _ in 0..1000 {
        // Ensure sqrt(a^2) = a or -a
        let a = Fq::rand(&mut rng);
        let mut nega = a;
        nega.negate();
        let mut b = a;
        b.square();

        let b = b.sqrt().unwrap();

        assert!(a == b || nega == b);
    }

    for _ in 0..1000 {
        // Ensure sqrt(a)^2 = a for random a
        let a = Fq::rand(&mut rng);

        if let Some(mut tmp) = a.sqrt() {
            tmp.square();

            assert_eq!(a, tmp);
        }
    }
}

#[test]
fn test_fq_from_into_repr() {
    // q should not be in the field
    assert!(Fq::from_repr(Fq::char()).is_err());

    // Zero should be in the field.
    assert!(Fq::from_repr(FqRepr::from(0)).unwrap().is_zero());

    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    for _ in 0..1000 {
        // Try to turn Fq elements into representations and back again, and compare.
        let a = Fq::rand(&mut rng);
        let a_repr = a.into_repr();
        let b_repr = FqRepr::from(a);
        assert_eq!(a_repr, b_repr);
        let a_again = Fq::from_repr(a_repr).unwrap();

        assert_eq!(a, a_again);
    }
}

#[test]
fn test_fq_num_bits() {
    println!("NUM_BITS, {}", Fq::NUM_BITS);
    println!("CAPACITY, {}", Fq::CAPACITY);
}

#[test]
fn fq_field_tests() {
    crate::tests::field::random_field_tests::<Fq>();
    crate::tests::field::random_sqrt_tests::<Fq>();
    crate::tests::field::random_frobenius_tests::<Fq, _>(Fq::char(), 13);
    crate::tests::field::from_str_tests::<Fq>();
}

#[test]
fn test_fq_ordering() {
    // FqRepr's ordering is well-tested, but we still need to make sure the Fq
    // elements aren't being compared in Montgomery form.
    for i in 0..100 {
        assert!(
            Fq::from_repr(FqRepr::from(i + 1)).unwrap() > Fq::from_repr(FqRepr::from(i)).unwrap()
        );
    }
}

#[test]
fn fq_repr_tests() {
    crate::tests::repr::random_repr_tests::<FqRepr>();
}

#[test]
fn test_fq_legendre() {
    use ff::LegendreSymbol::*;
    use ff::SqrtField;

    assert_eq!(QuadraticResidue, Fq::one().legendre());
    assert_eq!(Zero, Fq::zero().legendre());

    assert_eq!(
        QuadraticNonResidue,
        Fq::from_repr(FqRepr::from(11)).unwrap().legendre()
    );
    assert_eq!(
        QuadraticResidue,
        Fq::from_repr(FqRepr::from(4)).unwrap().legendre()
    );
}
