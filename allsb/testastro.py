

import logging
import json
from email.generator import BytesGenerator
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from email.encoders import encode_noop
from io import BytesIO
import os.path
import time

import astropy.io.fits as fits

import requests
import http.client
import allsb.utils as U

_logger = logging.getLogger(__name__)


AN_API_KEY = 'gmbziuhfdkfnpysg'
URL_API = "http://nova.astrometry.net/api"
URL_BASE = "http://nova.astrometry.net"


def wcs_calibrate_astrometry_net(filename, center_ij):
    # _logger.debug('shape is {shape}'.format(**datafile))
    _logger.debug('calibrate in astrometry.net')

    session = requests.Session()
    session_key = astro_session_key(session)

    upload_args = {
        'scale_upper': 90.0, 'publicly_visible': 'y',
        'allow_modifications': 'd', 'scale_type': 'ul',
        'allow_commercial_use': 'd', 'scale_lower': 25.0
    }
    upload_args['session'] = session_key

    response = astro_upload_file(session, filename, center_ij, upload_args=upload_args)

    url_sub = URL_API + "/submissions/{subid}"
    url_jobs = URL_API + "/jobs/{}"

    url = url_sub.format(**response.json())

    # Loop to check the status of the processing
    while True:
        try:
            response = session.get(url)
            mm = response.json()
            _logger.debug("response.json() %s", mm)
            time.sleep(10)
            _logger.debug('jobs: %s', mm['jobs'])
            if mm['jobs']:
                for jid in mm['jobs']:
                    if jid is not None:
                        response_job = session.get(url_jobs.format(jid))
                        _logger.debug('job %s %s', jid, response_job.text)
            _logger.debug('job_calibrations: %s', mm['job_calibrations'])
            if mm['job_calibrations']:
                break
        except http.client.RemoteDisconnected as error:
            _logger.error("%s", error)

    astro_download_files(session, jid, filename)


def astro_download_files(session, jobid, filename):


    url_results = URL_API + "/jobs/{}/calibration/"
    url_newfits = URL_BASE + "/new_fits_file/{}"

    url_corrfile = URL_BASE + "/corr_file/{}"

    response = session.get(url_results.format(jobid))
    nfile = U.insert_adj_filename(filename, 'results', ext='.json')
    _logger.debug('save astrometry.net results in %s', nfile)
    with open(nfile, 'wb') as fd:
        fd.write(response.content)

    urls = [url_newfits.format(jobid), url_corrfile.format(jobid)]
    ajs_s = ['newfits', 'corrfile']

    for url, adj in zip(urls, ajs_s):
        response = session.get(url)
        nfile = U.insert_adj_filename(filename, adj)
        _logger.debug('save in %s', nfile)
        with open(nfile, 'wb') as fd:
            fd.write(response.content)


def astro_session_key(session):

    rjson = json.dumps({"apikey": AN_API_KEY})
    url_login = URL_API + '/login'
    result = session.post(url_login, data={'request-json': rjson})
    _logger.debug('login %s', result.json())
    sess = result.json().get('session')
    _logger.debug('Got session: %s', sess)
    if not sess:
        raise ValueError('No "session" in response')
    return sess


def cut_center(filename, center, hsize=500):
    with fits.open(filename) as hdulist:
        data = hdulist[0].data
        _logger.debug('coordinates %s %s %s %s', center[0] - hsize, center[0] + hsize, center[1] - hsize, center[1] + hsize)

        newdata = data[center[0] - hsize:center[0] + hsize,
                  center[1] - hsize:center[1] + hsize]
        newhdu = fits.PrimaryHDU(newdata, hdulist[0].header)
        return newhdu


def calc_filename_center(filename):
    return U.insert_adj_filename(filename, 'center')


def astro_upload_file(session, filename, center_ij, hsize=500, upload_args=None):

    if upload_args is None:
        upload_args = {}

    _logger.debug('computing cut image around (i,j)=%s', center_ij)
    cut_hdu = cut_center(filename, center_ij, hsize)
    cut_filename = calc_filename_center(filename)
    cut_hdu.writeto(cut_filename, overwrite=True)

    try:
        with open(cut_filename, 'rb') as f:
            file_args = (cut_filename, f.read())
    except IOError as error:
        _logger.error("%s", error)
        raise

    headers, data = data_encode(upload_args, file_args)
    url_upload = URL_API + '/upload'
    response = session.post(url_upload, headers=headers, data=data)
    return response


class MyBytesGenerator(BytesGenerator):
    def __init__(self, fp, root=True):
        super().__init__(fp, mangle_from_=False,maxheaderlen=0)
        self.root = root

    def _write_headers(self, msg):
        # We don't want to write the top-level headers;
        # they go into Request(headers) instead.
        if self.root:
            return
        # We need to use \r\n line-terminator, but Generator
        # doesn't provide the flexibility to override, so we
        # have to copy-n-paste-n-modify.
        for h, v in msg.items():
            self._fp.write(('%s: %s\r\n' % (h,v)).encode())
        # A blank line always separates headers from body
        self._fp.write('\r\n'.encode())

    # The _write_multipart method calls "clone" for the
    # subparts.  We hijack that, setting root=False
    def clone(self, fp):
        return MyBytesGenerator(fp, root=False)


def data_encode(args, file_args):

    # Make a custom generator to format it the way we need.
    json_args = json.dumps(args)
    # If we're sending a file, format a multipart/form-data

    m1 = MIMEBase('text', 'plain')
    m1.add_header('Content-disposition',
                  'form-data; name="request-json"')
    m1.set_payload(json_args)
    m2 = MIMEApplication(file_args[1], 'octet-stream', encode_noop)
    m2.add_header('Content-disposition',
                  'form-data; name="file"; filename="%s"' % file_args[0])
    mp = MIMEMultipart('form-data', None, [m1, m2])

    fp = BytesIO()
    g = MyBytesGenerator(fp)
    g.flatten(mp)
    data = fp.getvalue()
    headers = {'Content-type': mp.get('Content-type')}

    return headers, data
