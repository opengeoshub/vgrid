import os
import zipfile
import tarfile
import json
import requests
from urllib.parse import urlparse
# Reference: https://github.com/opengeos/segment-geospatial/blob/main/samgeo/common.py


def github_raw_url(url):
    """Get the raw URL for a GitHub file.

    Args:
        url (str): The GitHub URL.
    Returns:
        str: The raw URL.
    """
    if isinstance(url, str) and url.startswith("https://github.com/") and "blob" in url:
        url = url.replace("github.com", "raw.githubusercontent.com").replace(
            "blob/", ""
        )
    return url

def install_package(package):
    """Install a Python package.

    Args:
        package (str | list): The package name or a GitHub URL or a list of package names or GitHub URLs.
    """
    import subprocess

    if isinstance(package, str):
        packages = [package]
    elif isinstance(package, list):
        packages = package

    for package in packages:
        if package.startswith("https"):
            package = f"git+{package}"

        # Execute pip install command and show output in real-time
        command = f"pip install {package}"
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)

        # Print output in real-time
        while True:
            output = process.stdout.readline()
            if output == b"" and process.poll() is not None:
                break
            if output:
                print(output.decode("utf-8").strip())

        # Wait for process to complete
        process.wait()

def extract_archive(archive, outdir=None, **kwargs) -> None:
    """
    Extracts a multipart archive.

    This function uses the patoolib library to extract a multipart archive.
    If the patoolib library is not installed, it attempts to install it.
    If the archive does not end with ".zip", it appends ".zip" to the archive name.
    If the extraction fails (for example, if the files already exist), it skips the extraction.

    Args:
        archive (str): The path to the archive file.
        outdir (str): The directory where the archive should be extracted.
        **kwargs: Arbitrary keyword arguments for the patoolib.extract_archive function.

    Returns:
        None

    Raises:
        Exception: An exception is raised if the extraction fails for reasons other than the files already existing.

    Example:

        files = ["sam_hq_vit_tiny.zip", "sam_hq_vit_tiny.z01", "sam_hq_vit_tiny.z02", "sam_hq_vit_tiny.z03"]
        base_url = "https://github.com/opengeos/datasets/releases/download/models/"
        urls = [base_url + f for f in files]
        leafmap.download_files(urls, out_dir="models", multi_part=True)

    """
    try:
        import patoolib
    except ImportError:
        install_package("patool")
        import patoolib

    if not archive.endswith(".zip"):
        archive = archive + ".zip"

    if outdir is None:
        outdir = os.path.dirname(archive)

    try:
        patoolib.extract_archive(archive, outdir=outdir, **kwargs)
    except Exception as e:
        print("The unzipped files might already exist. Skipping extraction.")
        return


def download_file(
    url=None,
    output=None,
    quiet=False,
    proxy=None,
    speed=None,
    use_cookies=True,
    verify=True,
    id=None,
    fuzzy=False,
    resume=False,
    unzip=True,
    overwrite=False,
    subfolder=False,
):
    """Download a file from URL, including Google Drive shared URL.

    Args:
        url (str, optional): Google Drive URL is also supported. Defaults to None.
        output (str, optional): Output filename. Default is basename of URL.
        quiet (bool, optional): Suppress terminal output. Default is False.
        proxy (str, optional): Proxy. Defaults to None.
        speed (float, optional): Download byte size per second (e.g., 256KB/s = 256 * 1024). Defaults to None.
        use_cookies (bool, optional): Flag to use cookies. Defaults to True.
        verify (bool | str, optional): Either a bool, in which case it controls whether the server's TLS certificate is verified, or a string,
            in which case it must be a path to a CA bundle to use. Default is True.. Defaults to True.
        id (str, optional): Google Drive's file ID. Defaults to None.
        fuzzy (bool, optional): Fuzzy extraction of Google Drive's file Id. Defaults to False.
        resume (bool, optional): Resume the download from existing tmp file if possible. Defaults to False.
        unzip (bool, optional): Unzip the file. Defaults to True.
        overwrite (bool, optional): Overwrite the file if it already exists. Defaults to False.
        subfolder (bool, optional): Create a subfolder with the same name as the file. Defaults to False.

    Returns:
        str: The output file path.
    """
    try:
        import gdown
    except ImportError:
        print(
            "The gdown package is required for this function. Use `pip install gdown` to install it."
        )
        return

    if output is None:
        if isinstance(url, str) and url.startswith("http"):
            output = os.path.basename(url)

    out_dir = os.path.abspath(os.path.dirname(output))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if isinstance(url, str):
        if os.path.exists(os.path.abspath(output)) and (not overwrite):
            print(
                f"{output} already exists. Skip downloading. Set overwrite=True to overwrite."
            )
            return os.path.abspath(output)
        else:
            url = github_raw_url(url)

    if "https://drive.google.com/file/d/" in url:
        fuzzy = True

    output = gdown.download(
        url, output, quiet, proxy, speed, use_cookies, verify, id, fuzzy, resume
    )

    if unzip:
        if output.endswith(".zip"):
            with zipfile.ZipFile(output, "r") as zip_ref:
                if not quiet:
                    print("Extracting files...")
                if subfolder:
                    basename = os.path.splitext(os.path.basename(output))[0]

                    output = os.path.join(out_dir, basename)
                    if not os.path.exists(output):
                        os.makedirs(output)
                    zip_ref.extractall(output)
                else:
                    zip_ref.extractall(os.path.dirname(output))
        elif output.endswith(".tar.gz") or output.endswith(".tar"):
            if output.endswith(".tar.gz"):
                mode = "r:gz"
            else:
                mode = "r"

            with tarfile.open(output, mode) as tar_ref:
                if not quiet:
                    print("Extracting files...")
                if subfolder:
                    basename = os.path.splitext(os.path.basename(output))[0]
                    output = os.path.join(out_dir, basename)
                    if not os.path.exists(output):
                        os.makedirs(output)
                    tar_ref.extractall(output)
                else:
                    tar_ref.extractall(os.path.dirname(output))

    return os.path.abspath(output)


def download_files(
    urls,
    out_dir=None,
    filenames=None,
    quiet=False,
    proxy=None,
    speed=None,
    use_cookies=True,
    verify=True,
    id=None,
    fuzzy=False,
    resume=False,
    unzip=True,
    overwrite=False,
    subfolder=False,
    multi_part=False,
):
    """Download files from URLs, including Google Drive shared URL.

    Args:
        urls (list): The list of urls to download. Google Drive URL is also supported.
        out_dir (str, optional): The output directory. Defaults to None.
        filenames (list, optional): Output filename. Default is basename of URL.
        quiet (bool, optional): Suppress terminal output. Default is False.
        proxy (str, optional): Proxy. Defaults to None.
        speed (float, optional): Download byte size per second (e.g., 256KB/s = 256 * 1024). Defaults to None.
        use_cookies (bool, optional): Flag to use cookies. Defaults to True.
        verify (bool | str, optional): Either a bool, in which case it controls whether the server's TLS certificate is verified, or a string, in which case it must be a path to a CA bundle to use. Default is True.. Defaults to True.
        id (str, optional): Google Drive's file ID. Defaults to None.
        fuzzy (bool, optional): Fuzzy extraction of Google Drive's file Id. Defaults to False.
        resume (bool, optional): Resume the download from existing tmp file if possible. Defaults to False.
        unzip (bool, optional): Unzip the file. Defaults to True.
        overwrite (bool, optional): Overwrite the file if it already exists. Defaults to False.
        subfolder (bool, optional): Create a subfolder with the same name as the file. Defaults to False.
        multi_part (bool, optional): If the file is a multi-part file. Defaults to False.

    Examples:

        files = ["sam_hq_vit_tiny.zip", "sam_hq_vit_tiny.z01", "sam_hq_vit_tiny.z02", "sam_hq_vit_tiny.z03"]
        base_url = "https://github.com/opengeos/datasets/releases/download/models/"
        urls = [base_url + f for f in files]
        leafmap.download_files(urls, out_dir="models", multi_part=True)
    """

    if out_dir is None:
        out_dir = os.getcwd()

    if filenames is None:
        filenames = [None] * len(urls)

    filepaths = []
    for url, output in zip(urls, filenames):
        if output is None:
            filename = os.path.join(out_dir, os.path.basename(url))
        else:
            filename = os.path.join(out_dir, output)

        filepaths.append(filename)
        if multi_part:
            unzip = False

        download_file(
            url,
            filename,
            quiet,
            proxy,
            speed,
            use_cookies,
            verify,
            id,
            fuzzy,
            resume,
            unzip,
            overwrite,
            subfolder,
        )

    if multi_part:
        archive = os.path.splitext(filename)[0] + ".zip"
        out_dir = os.path.dirname(filename)
        extract_archive(archive, out_dir)

        for file in filepaths:
            os.remove(file)


def download_folder(
    url=None,
    id=None,
    output=None,
    quiet=False,
    proxy=None,
    speed=None,
    use_cookies=True,
    remaining_ok=False,
):
    """Downloads the entire folder from URL.

    Args:
        url (str, optional): URL of the Google Drive folder. Must be of the format 'https://drive.google.com/drive/folders/{url}'. Defaults to None.
        id (str, optional): Google Drive's folder ID. Defaults to None.
        output (str, optional):  String containing the path of the output folder. Defaults to current working directory.
        quiet (bool, optional): Suppress terminal output. Defaults to False.
        proxy (str, optional): Proxy. Defaults to None.
        speed (float, optional): Download byte size per second (e.g., 256KB/s = 256 * 1024). Defaults to None.
        use_cookies (bool, optional): Flag to use cookies. Defaults to True.
        resume (bool, optional): Resume the download from existing tmp file if possible. Defaults to False.

    Returns:
        list: List of files downloaded, or None if failed.
    """

    try:
        import gdown
    except ImportError:
        print(
            "The gdown package is required for this function. Use `pip install gdown` to install it."
        )
        return

    files = gdown.download_folder(
        url, id, output, quiet, proxy, speed, use_cookies, remaining_ok
    )
    return files