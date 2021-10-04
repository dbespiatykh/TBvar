# %%
import os
import io
import pandas as pd
import numpy as np
import vcf
from gzip import open as gzopen

# %%
def vcf_to_dataframe(vcf_file):

    ext = os.path.splitext(vcf_file)[1]
    if ext == '.gz':
        file = gzopen(vcf_file, "rt")
    else:
        file = open(vcf_file)
    vcf_reader = vcf.Reader(file, 'r')
    res = []
    cols = ['sample', 'REF', 'ALT', 'pos']

    for rec in vcf_reader:
        x = [rec.end]
        for sample in rec.samples:
            if sample.gt_bases == None:
                # no call
                row = [sample.sample, rec.REF, sample.gt_bases]
            elif rec.REF != sample.gt_bases:
                cdata = sample.data
                row = [sample.sample, rec.REF, sample.gt_bases] + x
            else:
                # call is REF
                cdata = sample.data
                row = [sample.sample, rec.REF, sample.gt_bases] + x

            res.append(row)
    res = pd.DataFrame(res, columns=cols)
    res = res[~res.pos.isnull()]
    return res


# %%
df = vcf_to_dataframe(snakemake.input[0])

# %%
df['pos'] = df['pos'].astype(str).replace('\.0', '', regex=True)

# %%
conditions = [
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('615938')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('4404247')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('3021283')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('529363')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('3216553')),
    (df['REF'].eq('T') & df['ALT'].eq('C') & df['pos'].eq('1750465')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('2622402')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1491275')),
    (df['REF'].eq('A') & df['ALT'].eq('G') & df['pos'].eq('1136017')),
    (df['REF'].eq('C') & df['ALT'].eq('A') & df['pos'].eq('3479545')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('590595')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1273058')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('1923530')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('782943')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('1259841')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('2378976')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('3296809')),
    (df['REF'].eq('G') & df['ALT'].eq('C') & df['pos'].eq('1317597')),
    (df['REF'].eq('T') & df['ALT'].eq('C') & df['pos'].eq('2070728')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('962445')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('3470377')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('528781')),
    (df['REF'].eq('A') & df['ALT'].eq('C') & df['pos'].eq('2737201')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('2763624')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('1245275')),
    (df['REF'].eq('T') & df['ALT'].eq('C') & df['pos'].eq('61842')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('497491')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('282892')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('1881090')),
    (df['REF'].eq('A') & df['ALT'].eq('G') & df['pos'].eq('648465')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('2505085')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('797736')),
    (df['REF'].eq('C') & df['ALT'].eq('G') & df['pos'].eq('2078246')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('4248115')),
    (df['REF'].eq('A') & df['ALT'].eq('G') & df['pos'].eq('1947282')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('3836274')),
    (df['REF'].eq('A') & df['ALT'].eq('G') & df['pos'].eq('1692069')),
    (df['REF'].eq('G') & df['ALT'].eq('T') & df['pos'].eq('346693')),
    (df['REF'].eq('C') & df['ALT'].eq('A') & df['pos'].eq('3273107')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('342873')),
    (df['REF'].eq('C') & df['ALT'].eq('G') & df['pos'].eq('958362')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1084911')),
    (df['REF'].eq('G') & df['ALT'].eq('T') & df['pos'].eq('1591545')),
    (df['REF'].eq('G') & df['ALT'].eq('C') & df['pos'].eq('3722702')),
    (df['REF'].eq('C') & df['ALT'].eq('G') & df['pos'].eq('1237818')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('2874344')),
    (df['REF'].eq('G') & df['ALT'].eq('C') & df['pos'].eq('17842')),
    (df['REF'].eq('T') & df['ALT'].eq('C') & df['pos'].eq('931123')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('2825466')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('62657')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('514245')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('1850119')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1006080')),
    (df['REF'].eq('T') & df['ALT'].eq('G') & df['pos'].eq('541048')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1109535')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('4229087')),
    (df['REF'].eq('C') & df['ALT'].eq('A') & df['pos'].eq('286300')),
    (df['REF'].eq('A') & df['ALT'].eq('G') & df['pos'].eq('891756')),
    (df['REF'].eq('A') & df['ALT'].eq('G') & df['pos'].eq('3147742')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('107794')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('342340')),
    (df['REF'].eq('C') & df['ALT'].eq('G') & df['pos'].eq('2488724')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('1564799')),
    (df['REF'].eq('G') & df['ALT'].eq('C') & df['pos'].eq('58786')),
    (df['REF'].eq('G') & df['ALT'].eq('C') & df['pos'].eq('2411730')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('1466779')),
    (df['REF'].eq('A') & df['ALT'].eq('C') & df['pos'].eq('783601')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('870187')),
    (df['REF'].eq('C') & df['ALT'].eq('A') & df['pos'].eq('1487796')),
    (df['REF'].eq('T') & df['ALT'].eq('C') & df['pos'].eq('353766')),
    (df['REF'].eq('T') & df['ALT'].eq('C') & df['pos'].eq('1455780')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('611463')),
    (df['REF'].eq('C') & df['ALT'].eq('G') & df['pos'].eq('764995')),
    (df['REF'].eq('C') & df['ALT'].eq('A') & df['pos'].eq('1452071')),
    (df['REF'].eq('C') & df['ALT'].eq('A') & df['pos'].eq('615614')),
    (df['REF'].eq('T') & df['ALT'].eq('C') & df['pos'].eq('825585')),
    (df['REF'].eq('T') & df['ALT'].eq('C') & df['pos'].eq('1647807')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('4316114')),
    (df['REF'].eq('G') & df['ALT'].eq('C') & df['pos'].eq('3414791')),
    (df['REF'].eq('C') & df['ALT'].eq('G') & df['pos'].eq('3388166')),
    (df['REF'].eq('G') & df['ALT'].eq('C') & df['pos'].eq('784581')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('403364')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('2077253')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('3977226')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1297327')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('4398141')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1274335')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('1132368')),
    (df['REF'].eq('G') & df['ALT'].eq('T') & df['pos'].eq('784440')),
    (df['REF'].eq('C') & df['ALT'].eq('A') & df['pos'].eq('1502120')),
    (df['REF'].eq('T') & df['ALT'].eq('C') & df['pos'].eq('225495')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('4307886')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('4151558')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('2905505')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('355181')),
    (df['REF'].eq('C') & df['ALT'].eq('G') & df['pos'].eq('15036')),
    (df['REF'].eq('G') & df['ALT'].eq('C') & df['pos'].eq('2694560')),
    (df['REF'].eq('C') & df['ALT'].eq('G') & df['pos'].eq('342201')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('4246508')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('985287')),
    (df['REF'].eq('G') & df['ALT'].eq('T') & df['pos'].eq('1719757')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('620029')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('3466426')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('18091')),
    (df['REF'].eq('G') & df['ALT'].eq('C') & df['pos'].eq('4260268')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('874787')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('4406749')),
    (df['REF'].eq('G') & df['ALT'].eq('C') & df['pos'].eq('1501468')),
    (df['REF'].eq('C') & df['ALT'].eq('G') & df['pos'].eq('1098698')),
    (df['REF'].eq('G') & df['ALT'].eq('C') & df['pos'].eq('4125058')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('4260742')),
    (df['REF'].eq('C') & df['ALT'].eq('G') & df['pos'].eq('3570528')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('896119')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('2875883')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('734562')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('4236903')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('17665')),
    (df['REF'].eq('C') & df['ALT'].eq('G') & df['pos'].eq('4249732')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('716918')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('3836739')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1130526')),
    (df['REF'].eq('G') & df['ALT'].eq('C') & df['pos'].eq('2914906')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('2417281')),
    (df['REF'].eq('C') & df['ALT'].eq('G') & df['pos'].eq('616408')),
    (df['REF'].eq('G') & df['ALT'].eq('T') & df['pos'].eq('1759252')),
    (df['REF'].eq('A') & df['ALT'].eq('G') & df['pos'].eq('420008')),
    (df['REF'].eq('C') & df['ALT'].eq('G') & df['pos'].eq('119600')),
    (df['REF'].eq('C') & df['ALT'].eq('A') & df['pos'].eq('1799921')),
    (df['REF'].eq('C') & df['ALT'].eq('G') & df['pos'].eq('1816587')),
    (df['REF'].eq('G') & df['ALT'].eq('T') & df['pos'].eq('982363')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1137518')),
    (df['REF'].eq('G') & df['ALT'].eq('T') & df['pos'].eq('221190')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('2750052')),
    (df['REF'].eq('A') & df['ALT'].eq('G') & df['pos'].eq('2831482')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('1882180')),
    (df['REF'].eq('A') & df['ALT'].eq('G') & df['pos'].eq('62768')),
    (df['REF'].eq('T') & df['ALT'].eq('C') & df['pos'].eq('904090')),
    (df['REF'].eq('G') & df['ALT'].eq('T') & df['pos'].eq('44812')),
    (df['REF'].eq('A') & df['ALT'].eq('C') & df['pos'].eq('2955957')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('1477596')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('1881090')),
    (df['REF'].eq('G') & df['ALT'].eq('T') & df['pos'].eq('346693')),
    (df['REF'].eq('T') & df['ALT'].eq('C') & df['pos'].eq('886115')),
    (df['REF'].eq('C') & df['ALT'].eq('A') & df['pos'].eq('1477522')),
    (df['REF'].eq('A') & df['ALT'].eq('G') & df['pos'].eq('782634')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('2532616')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1059643')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('2149762')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1364706')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('1518365')),
    (df['REF'].eq('G') & df['ALT'].eq('T') & df['pos'].eq('4086')),
    (df['REF'].eq('A') & df['ALT'].eq('G') & df['pos'].eq('1844339')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1329881')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1772800')),
    (df['REF'].eq('A') & df['ALT'].eq('C') & df['pos'].eq('2323004')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1479479')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('3459402')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('870862')),
    (df['REF'].eq('A') & df['ALT'].eq('G') & df['pos'].eq('1606119')),
    (df['REF'].eq('A') & df['ALT'].eq('G') & df['pos'].eq('162878')),
    (df['REF'].eq('C') & df['ALT'].eq('G') & df['pos'].eq('3680254')),
    (df['REF'].eq('T') & df['ALT'].eq('C') & df['pos'].eq('1182150')),
    (df['REF'].eq('C') & df['ALT'].eq('A') & df['pos'].eq('202317')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('3366367')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('1331291')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1851172')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('214609')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('2411421')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('3278762')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('566189')),
    (df['REF'].eq('A') & df['ALT'].eq('G') & df['pos'].eq('1194516')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('2810437')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('581990')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1207038')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('152602')),
    (df['REF'].eq('A') & df['ALT'].eq('C') & df['pos'].eq('2574591')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('260751')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('63781')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('813155')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('1509067')),
    (df['REF'].eq('C') & df['ALT'].eq('A') & df['pos'].eq('199859')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('559350')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('1249771')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('295040')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('793892')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('403481')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('101902')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('264535')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('658500')),
    (df['REF'].eq('G') & df['ALT'].eq('T') & df['pos'].eq('285096')),
    (df['REF'].eq('G') & df['ALT'].eq('T') & df['pos'].eq('1158248')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('236217')),
    (df['REF'].eq('G') & df['ALT'].eq('A') & df['pos'].eq('657099')),
    (df['REF'].eq('C') & df['ALT'].eq('T') & df['pos'].eq('2738352')),
]

values = [
    'lineage1::615938G->A',
    'lineage1.1::4404247G->A',
    'lineage1.1.1::3021283G->A',
    'lineage1.1.1 (Clark et al.)::529363C->T',
    'lineage1.1.1.1::3216553G->A',
    'lineage1.1.1.1 (Clark et al.)::1750465T->C',
    'lineage1.1.2::2622402G->A',
    'lineage1.1.3::1491275G->A',
    'lineage1.2 (Clark et al.)::1136017A->G',
    'lineage1.2.1::3479545C->A',
    'lineage1.2.1 (Clark et al.)::590595G->A',
    'lineage1.2.1.1::1273058G->A',
    'lineage1.2.1.1::1923530C->T',
    'lineage1.2.1.1::782943G->A',
    'lineage1.2.1.2::1259841C->T',
    'lineage1.2.1.2::2378976C->T',
    'lineage1.2.1.2::3296809G->A',
    'lineage1.2.1.3::1317597G->C',
    'lineage1.2.1.3::2070728T->C',
    'lineage1.2.1.3::962445G->A',
    'lineage1.2.2::3470377C->T',
    'lineage1.2.2 (Clark et al.)::528781G->A',
    'lineage1.2.2.1 (Clark et al.)::2737201A->C',
    'lineage1.3 (Clark et al.)::2763624G->A',
    'lineage1.3.1 (Clark et al.)::1245275C->T',
    'lineage1.3.2 (Clark et al.)::61842T->C',
    'lineage2::497491G->A',
    'lineage2 (Clark et al.)::282892C->T',
    'lineage2.1::1881090C->T',
    'lineage2.1 (Clark et al.)::648465A->G',
    'lineage2.2::2505085G->A',
    'lineage2.2.1::797736C->T',
    'lineage2.2.1 (Clark et al.)::2078246C->G',
    'lineage2.2.1.1::4248115C->T',
    'lineage2.2.1.1 (Clark et al.)::1947282A->G',
    'lineage2.2.1.2::3836274G->A',
    'lineage2.2.1.2 (Clark et al.)::1692069A->G',
    'lineage2.2.2::346693G->T',
    'lineage3::3273107C->A',
    'lineage3 (Clark et al.)::342873C->T',
    'lineage3.1 (Clark et al.)::958362C->G',
    'lineage3.1.1::1084911G->A',
    'lineage3.1.1 (Clark et al.)::1591545G->T',
    'lineage3.1.2::3722702G->C',
    'lineage3.1.2.1::1237818C->G',
    'lineage3.1.2.2::2874344G->A',
    'lineage3.2 (Clark et al.)::17842G->C',
    'lineage4::931123T->C',
    'lineage4 (Clark et al.)::2825466G->A',
    'lineage4.1::62657G->A',
    'lineage4.1.1::514245C->T',
    'lineage4.1.1.1::1850119C->T',
    'lineage4.1.1.1 (Clark et al.)::1006080G->A',
    'lineage4.1.1.2::541048T->G',
    'lineage4.1.1.2 (Clark et al.)::1109535G->A',
    'lineage4.1.1.3::4229087C->T',
    'lineage4.1.1.3.1 (Clark et al.)::286300C->A',
    'lineage4.1.2::891756A->G',
    'lineage4.1.2 (Clark et al.)::3147742A->G',
    'lineage4.1.2.1::107794C->T',
    'lineage4.1.2.1 (Clark et al.)::342340C->T',
    'lineage4.1.2.1.1 (Clark et al.)::2488724C->G',
    'lineage4.1.3 (Clark et al.)::1564799C->T',
    'lineage4.1.4 (Clark et al.)::58786G->C',
    'lineage4.2::2411730G->C',
    'lineage4.2 (Clark et al.)::1466779C->T',
    'lineage4.2.1::783601A->C',
    'lineage4.2.1.1 (Clark et al.)::870187C->T',
    'lineage4.2.2::1487796C->A',
    'lineage4.2.2 (Clark et al.)::353766T->C',
    'lineage4.2.2.1::1455780T->C',
    'lineage4.2.2.2 (Clark et al.)::611463G->A',
    'lineage4.3::764995C->G',
    'lineage4.3 (Clark et al.)::1452071C->A',
    'lineage4.3.1::615614C->A',
    'lineage4.3.1 (Clark et al.)::825585T->C',
    'lineage4.3.1.1 (Clark et al.)::1647807T->C',
    'lineage4.3.2::4316114G->A',
    'lineage4.3.2 (Clark et al.)::3414791G->C',
    'lineage4.3.2.1::3388166C->G',
    'lineage4.3.2.1 (Clark et al.)::784581G->C',
    'lineage4.3.3::403364G->A',
    'lineage4.3.3 (Clark et al.)::2077253G->A',
    'lineage4.3.4::3977226G->A',
    'lineage4.3.4 (Clark et al.)::1297327G->A',
    'lineage4.3.4.1::4398141G->A',
    'lineage4.3.4.1 (Clark et al.)::1274335G->A',
    'lineage4.3.4.2::1132368C->T',
    'lineage4.3.4.2 (Clark et al.)::784440G->T',
    'lineage4.3.4.2.1::1502120C->A',
    'lineage4.3.4.2.1 (Clark et al.)::225495T->C',
    'lineage4.4::4307886G->A',
    'lineage4.4.1::4151558G->A',
    'lineage4.4.1 (Clark et al.)::2905505G->A',
    'lineage4.4.1.1::355181G->A',
    'lineage4.4.1.1.1 (Clark et al.)::15036C->G',
    'lineage4.4.1.2::2694560G->C',
    'lineage4.4.1.2 (Clark et al.)::342201C->G',
    'lineage4.4.2::4246508G->A',
    'lineage4.4.2 (Clark et al.)::985287G->A',
    'lineage4.5::1719757G->T',
    'lineage4.5 (Clark et al.)::620029C->T',
    'lineage4.6::3466426G->A',
    'lineage4.6 (Clark et al.)::18091G->A',
    'lineage4.6.1::4260268G->C',
    'lineage4.6.1.1::874787G->A',
    'lineage4.6.1.1 (Clark et al.)::4406749G->A',
    'lineage4.6.1.2::1501468G->C',
    'lineage4.6.1.2 (Clark et al.)::1098698C->G',
    'lineage4.6.2::4125058G->C',
    'lineage4.6.2 (Clark et al.)::4260742G->A',
    'lineage4.6.2.1::3570528C->G',
    'lineage4.6.2.1 (Clark et al.)::896119C->T',
    'lineage4.6.2.2::2875883C->T',
    'lineage4.6.3 (Clark et al.)::734562G->A',
    'lineage4.6.4 (Clark et al.)::4236903G->A',
    'lineage4.6.5 (Clark et al.)::17665G->A',
    'lineage4.7::4249732C->G',
    'lineage4.7 (Clark et al.)::716918G->A',
    'lineage4.8::3836739G->A',
    'lineage4.8 (Clark et al.)::1130526G->A',
    'lineage4.8.1 (Clark et al.)::2914906G->C',
    'lineage4.8.2 (Clark et al.)::2417281G->A',
    'lineage4.8.3 (Clark et al.)::616408C->G',
    'lineage4.9::1759252G->T',
    'lineage4.9 (Clark et al.)::420008A->G',
    'lineage4.9.1 (Clark et al.)::119600C->G',
    'lineage5::1799921C->A',
    'lineage6::1816587C->G',
    'lineage6 (Clark et al.)::982363G->T',
    'lineage7::1137518G->A',
    'lineage8 (Clark et al.)::221190G->T',
    'lineage9 (Clark et al.)::2750052G->A',
    'lineageBOV::2831482A->G',
    'lineageBOV_AFRI::1882180C->T',
    'M.bovis (Clark et al.)::62768A->G',
    'M.caprae (Clark et al.)::904090T->C',
    'M.orygis (Clark et al.)::44812G->T',
    'Beijing::2955957A->C',
    'modern_Beijing::1477596C->T',
    'proto-beijing::1881090C->T',
    'Asia_Ancestral_1::346693G->T',
    'Asia_Ancestral_2::886115T->C',
    'Asia_Ancestral_3::1477522C->A',
    'Asian_African_1::782634A->G',
    'Asian_African_2::2532616G->A',
    'Asian_African_3::1059643G->A',
    'CAO::2149762G->A',
    'Central_Asia::1364706G->A',
    'CladeA::1518365C->T',
    'Europe/Russia_W148_outbreak::4086G->T',
    'Pacific_RD150::1844339A->G',
    'lineage1.1.1.2::1329881G->A',
    'lineage1.1.1.2::1772800G->A',
    'lineage1.1.1.2::2323004A->C',
    'lineage1.1.1.3::1479479G->A',
    'lineage1.1.1.3::3459402C->T',
    'lineage1.1.1.3::870862C->T',
    'lineage1.1.1.4::1606119A->G',
    'lineage1.1.1.4::162878A->G',
    'lineage1.1.1.4::3680254C->G',
    'lineage1.1.1.5::1182150T->C',
    'lineage1.1.1.5::202317C->A',
    'lineage1.1.1.5::3366367C->T',
    'lineage1.1.1.6::1331291C->T',
    'lineage1.1.1.6::1851172G->A',
    'lineage1.1.1.6::214609C->T',
    'lineage1.1.1.7::2411421C->T',
    'lineage1.1.1.7::3278762C->T',
    'lineage1.1.1.7::566189G->A',
    'lineage1.1.1.8::1194516A->G',
    'lineage1.1.1.8::2810437G->A',
    'lineage1.1.1.8::581990C->T',
    'lineage1.1.1.9::1207038G->A',
    'lineage1.1.1.9::152602C->T',
    'lineage1.1.1.9::2574591A->C',
    'lineage1.1.2.1::260751C->T',
    'lineage1.1.2.1::63781G->A',
    'lineage1.1.2.1::813155G->A',
    'lineage1.1.2.2::1509067C->T',
    'lineage1.1.2.2::199859C->A',
    'lineage1.1.2.2::559350G->A',
    'lineage1.1.3.1::1249771G->A',
    'lineage1.1.3.1::295040G->A',
    'lineage1.1.3.1::793892G->A',
    'lineage1.1.3.1 (Clark et al.)::403481C->T',
    'lineage1.1.3.2::101902C->T',
    'lineage1.1.3.2::264535C->T',
    'lineage1.1.3.2::658500C->T',
    'lineage1.1.3.2 (Clark et al.)::285096G->T',
    'lineage1.1.3.3::1158248G->T',
    'lineage1.1.3.3::236217C->T',
    'lineage1.1.3.3::657099G->A',
    'lineage1.1.3.3 (Clark et al.)::2738352C->T',
]

# %%
df['lineage'] = np.select(conditions, values)

# %%
df.drop(['REF', 'ALT', 'pos'], axis=1, inplace=True)

# %%
df = (
    df.groupby(['sample'])
    .agg(lambda x: ','.join(set(x)))
    .reset_index()
    .reindex(columns=df.columns)
)

# %%
df = df.join(df['lineage'].str.get_dummies(sep=',').astype(int))

# %%
df = df.drop(['lineage', '0'], axis=1).replace({1: "+", 0: ""})

# %%
def lineage2_decision(df):
    if (
        'Beijing::2955957A->C' in df
        and 'proto-beijing::1881090C->T' in df
        and df['Beijing::2955957A->C'] == "+"
        and df['proto-beijing::1881090C->T'] == "+"
    ):
        return 'proto-beijing::1881090C->T'
    elif (
        'Beijing::2955957A->C' in df
        and 'Asia_Ancestral_1::346693G->T' in df
        and df['Beijing::2955957A->C'] == "+"
        and df['Asia_Ancestral_1::346693G->T'] == "+"
    ):
        return 'Asia_Ancestral_1::346693G->T'
    elif (
        'Beijing::2955957A->C' in df
        and 'Asia_Ancestral_2::886115T->C' in df
        and df['Beijing::2955957A->C'] == "+"
        and df['Asia_Ancestral_2::886115T->C'] == "+"
    ):
        return 'Asia_Ancestral_2::886115T->C'
    elif (
        'Beijing::2955957A->C' in df
        and 'Asia_Ancestral_3::1477522C->A' in df
        and df['Beijing::2955957A->C'] == "+"
        and df['Asia_Ancestral_3::1477522C->A'] == "+"
    ):
        return 'Asia_Ancestral_3::1477522C->A'
    elif (
        'Beijing::2955957A->C' in df
        and 'Asian_African_1::782634A->G' in df
        and df['Beijing::2955957A->C'] == "+"
        and df['Asian_African_1::782634A->G'] == "+"
    ):
        return 'Asian_African_1::782634A->G'
    elif (
        'Beijing::2955957A->C' in df
        and 'Asian_African_2::2532616G->A' in df
        and df['Beijing::2955957A->C'] == "+"
        and df['Asian_African_2::2532616G->A'] == "+"
    ):
        return 'Asian_African_2::2532616G->A'
    elif (
        'Beijing::2955957A->C' in df
        and 'Asian_African_3::1059643G->A' in df
        and df['Beijing::2955957A->C'] == "+"
        and df['Asian_African_3::1059643G->A'] == "+"
    ):
        return 'Asian_African_3::1059643G->A'
    elif (
        'Beijing::2955957A->C' in df
        and 'CAO::2149762G->A' in df
        and df['Beijing::2955957A->C'] == "+"
        and df['CAO::2149762G->A'] == "+"
    ):
        return 'CAO::2149762G->A'
    elif (
        'Beijing::2955957A->C' in df
        and 'Central_Asia::1364706G->A' in df
        and df['Beijing::2955957A->C'] == "+"
        and df['Central_Asia::1364706G->A'] == "+"
    ):
        return 'Central_Asia::1364706G->A'
    elif (
        'Beijing::2955957A->C' in df
        and 'CladeA::1518365C->T' in df
        and df['Beijing::2955957A->C'] == "+"
        and df['CladeA::1518365C->T'] == "+"
    ):
        return 'CladeA::1518365C->T'
    elif (
        'Beijing::2955957A->C' in df
        and 'Europe/Russia_W148_outbreak::4086G->T' in df
        and df['Beijing::2955957A->C'] == "+"
        and df['Europe/Russia_W148_outbreak::4086G->T'] == "+"
    ):
        return 'Europe/Russia_W148_outbreak::4086G->T'
    elif (
        'Beijing::2955957A->C' in df
        and 'Pacific_RD150::1844339A->G' in df
        and df['Beijing::2955957A->C'] == "+"
        and df['Pacific_RD150::1844339A->G'] == "+"
    ):
        return 'Pacific_RD150::1844339A->G'
    elif 'Beijing::2955957A->C' not in df or df['Beijing::2955957A->C'] != "+":
        return np.nan


# %%
df["L2"] = df.apply(lineage2_decision, axis=1)

# %%
df.to_csv(snakemake.output[0], sep='\t', index=False)
