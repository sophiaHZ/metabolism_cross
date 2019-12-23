#!/usr/bin/env python3
# -*-coding:utf8-*-
# author: Todzhu
# Date: 2019/4/1 11:29

import os,argparse
from bs4 import BeautifulSoup
import threading, queue, time, requests, re

head = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/68.0.3440.106 Safari/537.36'}
parser = argparse.ArgumentParser(prog='KEGG_color_map.py',description="Use KEGG API to create color pathway map, then download the color png file and html")
parser.add_argument('-i', action='store', dest='map_info', help='Input the map information file: five columns, 1.mapId 2.protein id 3.protein to KO id 4.color 5.protein ratio 6.protein regulated type')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parse = parser.parse_args()

pathwayMap = {}
ko2prot = {}
baseUrl = 'https://www.kegg.jp/kegg-bin/show_pathway?'
urlQueue = queue.Queue()

# 新建image文件夹，用于存放pathway map png文件
if not os.path.exists('./image'):
    os.mkdir('./image')

# 构造url地址
def makeUrl(pathinfo):
    with open(pathinfo,'r') as f:
        for line in f.readlines():
            line = line.strip().split('\t')
            #ko2prot[line[2]] = line[2]+' ('+line[1]+' | '+line[4]+' '+line[5]+')'
            ko2prot.setdefault(line[2],[]).append(line[2]+' ('+line[1]+' | '+line[4]+' '+line[5]+')')
            pathwayMap.setdefault(line[0],[]).append(line[2]+'%09'+line[3])
        for key,value in pathwayMap.items():
            url = baseUrl+key+'/'+'+'.join(value)
            urlQueue.put(url)

# 访问url地址，抓取图片及网页内容
def fetchUrl(urlQueue):
    while True:
        time.sleep(2)
        try:
            #不阻塞的读取队列数据
            url = urlQueue.get_nowait()
            i = urlQueue.qsize()
        except Exception as e:
            break

        print ('Current Thread Name %s, Url: %s ' % (threading.currentThread().name, url))

        try:
            response = requests.get(url,headers=head,timeout=100)
            responseCode = response.status_code
            time.sleep(1.5)			
        except Exception as e:
            continue

        if responseCode == 200:
            #抓取内容的数据处理可以放到这里
            soup = BeautifulSoup(response.content, 'html.parser')
            map_url = 'https://www.kegg.jp' + soup.find('img', attrs={'name': 'pathwayimage'}).get('src')

            #下载KEGG pathway map PNG
            with open(os.path.join('./image',map_url.split('/')[-1]),'wb') as p:
                p.write(requests.get(map_url).content)

            with open(url.split('/')[4].split('?')[1]+'.html', 'w') as html:
                # replace KO ID to protein information
                for i in soup.find_all('area'):
                    title = re.findall('[KCGD]\d{5}',i['title'])
                    node = []
                    for key, value in ko2prot.items():
                        if key in title:
                            node.append('\n'.join(list(set(value))))
                        if node:
                            i['title'] = '\n'.join(node)
                            #i['title'] = re.sub(key + r' \(.*?\)', '\n'.join(value), i['title'])
							#i['title'] = '\n'.join(list(set(ko2prot[key])))

                response = soup.prettify().replace('href="/', 'href="https://www.kegg.jp/').replace('src="/','src="https://www.kegg.jp/'). \
                    replace('<img border="0" name="pathwayimage" src="%s" usemap="#mapdata"/>' % map_url,
                            '<img border="0" name="pathwayimage" src="image/%s" usemap="#mapdata"/>' % map_url.split('/')[-1]). \
                    replace('/kegg/document/help_pathway.html', 'https://www.kegg.jp/kegg/document/help_pathway.html').\
                    replace('align="middle" alt="Help" border="0" name="help" onmousedown="btn(this,\'Hbd\')" onmouseout="btn(this,\'Hb\')" onmouseover="btn(this,\'Hbh\')" onmouseup="btn(this,\'Hb\')"','')
                
                html.write(response.replace('\u2013','-'))
            # 设置延时
            time.sleep(2.5)

        else:
            print('下载失败：'+url)
            break


if __name__ == '__main__':
    startTime = time.time()
    makeUrl(parse.map_info)
    threads = []
    # 可以调节线程数， 进而控制抓取速度
    threadNum = 10
    for i in range(0, threadNum):
        t = threading.Thread(target=fetchUrl, args=(urlQueue,))
        threads.append(t)
    for t in threads:
        t.start()
    for t in threads:
        #多线程多join的情况下，依次执行各线程的join方法, 这样可以确保主线程最后退出， 且各个线程间没有阻塞
        t.join()
    endTime = time.time()
    print ('Done, Time cost: %s ' %  (endTime - startTime))
