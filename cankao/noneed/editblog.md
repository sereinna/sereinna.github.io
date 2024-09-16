---
title: fuwari 博客修改教程
published: 2024-04-01
description: ""
image: ""
tags: [Fuwari]
category: 指南
draft: false
lang: ""
---

## 修改 astro.config.mjs

```diff
- site: "https://fuwari.vercel.app/",
+ site: "https://your.domain",
```

## 修改./src/config.ts

```diff
- title: 'Fuwari',
- subtitle: 'Demo Site',
- lang: 'en',
+ title: 'your title',
+ subtitle: 'your subtitle',
+ lang: 'zh_CN',

- enable: false,
- src: 'assets/images/demo-banner.png',
+ enable: true,
+ src: './assets/banner.jpg', //添加相对路径图片

- enable: false,
- text: '',
- url: ''
+ enable: true,
+ text: '空色天絵 / NEO TOKYO NOIR 01',
+ url: 'https://www.pixiv.net/artworks/111024784'

- name: 'Lorem Ipsum',
- bio: 'Lorem ipsum dolor sit amet, consectetur adipiscing elit.',
+ name: 'profile name',
+ bio: 'profile desc',

# 修改、添加联系方式
  {
      name: 'GitHub',
      icon: 'fa6-brands:github',
      icon: 'fa6-brands:github',
-     url: 'https://github.com/saicaca/fuwari',
+     url: 'https://github.com/yourname',
+ },
+ {
+     name: 'Bilibili',
+     icon: 'fa6-brands:bilibili',
+     url: 'https://bilibili.com',
+ },
```

## 创建自定义页面(友链为例)

### ./src/types/config.ts

```diff
 export enum LinkPreset {
   Home = 0,
   Archive = 1,
   About = 2,
+  Links = 3,
 }
```

### ./src/config.ts

```diff
 export const navBarConfig: NavBarConfig = {
   links: [
     LinkPreset.Home,
     LinkPreset.Archive,
+     LinkPreset.Links,
     LinkPreset.About,
     {
       name: 'GitHub',
       url: 'https://github.com/saicaca/fuwari',
       external: true,
     },
   ],
 }
```

### ./src/constants/link-presets.ts

```diff
     name: i18n(I18nKey.archive),
     url: '/archive/',
   },
+  [LinkPreset.Links]: {
+    name: i18n(I18nKey.links),
+    url: '/links/',
+  }
 }
```

### ./src/i18n/i18nKey.ts

```diff
  home = 'home',
  about = 'about',
  archive = 'archive',
+ links = 'links',
  search = 'search',
```

### ./src/i18n/languages/zh_CN.ts

```diff
  [Key.home]: '主页',
  [Key.about]: '关于',
  [Key.archive]: '归档',
+ [Key.links]: '友链',
  [Key.search]: '搜索',
```

### ./src/pages/links.astro

```js
---
import MainGridLayout from '../layouts/MainGridLayout.astro'

import { getEntry } from 'astro:content'
import Markdown from '@components/misc/Markdown.astro'
import I18nKey from '@i18n/i18nKey'
import { i18n } from '@i18n/translation'
// 导入links.md
const linksPost = await getEntry('spec', 'links')
const { Content } = await linksPost.render()
const items = [
  {
    title: 'Astro',
    imgurl: 'https://docs.astro.build/favicon.ico',
    desc: 'Astro',
    siteurl: 'https://astro.build/',
    tags: ['框架'],
  }
]
---
<MainGridLayout title={i18n(I18nKey.links)} description={i18n(I18nKey.links)}>
    <div class="flex w-full rounded-[var(--radius-large)] overflow-hidden relative min-h-32">
        <!-- 链接 -->
        <div class="card-base z-10 px-9 py-6 relative w-full ">
            <div class="grid grid-cols-1 sm:grid-cols-2 gap-x-6 gap-y-8 my-4">
                {items.map((item) => (
                    <div class="flex flex-nowrap items-stretch h-28 gap-4 rounded-[var(--radius-large)]">
                        <div class="w-28 h-28 flex-shrink-0 rounded-lg overflow-hidden bg-zinc-200 dark:bg-zinc-900">
                            <img src={item.imgurl} alt="urlimg" class="w-full h-full object-cover">
                        </div>
                        <div class="grow w-full">
                            <div class="font-bold transition text-lg text-neutral-900 dark:text-neutral-100 mb-1">{item.title}</div>
                            <div class="text-50 text-sm font-medium">{item.desc}</div>
                            <div class:list={["items-center", {"flex": true, "hidden md:flex" : false}]}>
                                <div class="flex flex-row flex-nowrap items-center">
                                    {(item.tags && item.tags.length > 0) && item.tags.map((tag,i) => (
                                    <div class:list={[{"hidden": i==0}, "mx-1.5 text-[var(--meta-divider)] text-sm" ]}>
                                        /
                                    </div>
                                    <span class="transition text-50 text-sm font-medium">
                                        {tag}
                                    </span>))}
                                    {!(item.tags && item.tags.length > 0) && <div class="transition text-50 text-sm font-medium">{i18n(I18nKey.noTags)}</div>}
                                </div>
                            </div>
                        </div>
                        <a href={item.siteurl} target="_blank" rel="noopener noreferrer"class="flex btn-regular w-[3.25rem] rounded-lg bg-[var(--enter-btn-bg)] hover:bg-[var(--enter-btn-bg-hover)] active:bg-[var(--enter-btn-bg-active)] active:scale-95">
                            <svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" aria-hidden="true" role="img" class="transition text-[var(--primary)] text-4xl mx-auto iconify iconify--material-symbols" width="1em" height="1em" viewBox="0 0 24 24">
                                <path fill="currentColor" d="M12.6 12L8.7 8.1q-.275-.275-.275-.7t.275-.7t.7-.275t.7.275l4.6 4.6q.15.15.213.325t.062.375t-.062.375t-.213.325l-4.6 4.6q-.275.275-.7.275t-.7-.275t-.275-.7t.275-.7z"></path>
                            </svg>
                        </a>
                    </div>
                ))}
            </div>
            <Markdown class="mt-2">
                <Content />
            </Markdown>
        </div>
    </div>
</MainGridLayout>
```

### ./src/content/spes/links.md

```markdown
# FriendLinks

欢迎各位
```
